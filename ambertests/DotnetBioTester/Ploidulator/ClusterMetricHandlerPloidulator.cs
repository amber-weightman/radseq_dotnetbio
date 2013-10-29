using Bio.Algorithms.Alignment;
using Bio.Algorithms.Metric;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.DataVisualization.Charting;
using System.Windows.Data;
using System.Windows.Media;
using System.Windows.Threading;

namespace Ploidulator
{
    /// <summary>
    /// An ClusterMetricHandler receives a set of SAMAlignedSequences (as a list or one by one)
    /// and passes them as "clusters" to a ClusterMetric, which performs a number of calculations
    /// on them to calculate the accuracy of the cluster.
    /// </summary>
    public class ClusterMetricHandlerPloidulator : IMetricHandler
    {
        #region Private Static Fields

        private static string OUTPUT_FOLDER = @"E:\Harvard\pl_output";
        private static int OUTPUT_QUEUE_SIZE = 7; // max number of sequences that can be stored in the
                                                     // output queue (to prevent too many being held in memory)
        #endregion

        #region Private Fields

        /// <summary>
        /// Name of file to which metric data will be written.
        /// </summary>
        private string fileName;

        /// <summary>
        /// The expected ploidy of the organism being sequenced, e.g. 2 (diploid), 3 (triploid)
        /// </summary>
        private int expectedPloidy;

        /// <summary>
        /// Whether the handler has been aborted
        /// </summary>
        private bool aborted = false;

        /// <summary>
        /// Whether the handler has done its completion cleanup 
        /// </summary>
        private bool isComplete = false;

        /// <summary>
        /// Id of the cluster to which current sequences belong
        /// </summary>
        private string currCluster;

        /// <summary>
        /// Number of clusters which have been received
        /// </summary>
        private int clusterCount = 0;

        /// <summary>
        /// Original header of input file
        /// </summary>
        private SAMAlignmentHeader header = null;

        /// <summary>
        /// New header for filtered output file
        /// </summary>
        private SAMAlignmentHeader newHeader = null;

        #region flags
        
        /// <summary>
        /// Assuming a single bam output file and a single thread writing to this file, this value is true
        /// if it is ok to write to the bam file
        /// </summary>
        private bool canWriteToBam = true;

        /// <summary>
        /// Whether the sequence and header bam output files have been merged
        /// </summary>
        private bool bamFilesMerged = false;
        
        /// <summary>
        /// This value is true if processing has finished for the current input file
        /// </summary>
        private bool finished = false;

        #endregion

        #region storage for sequences

        /// <summary>
        /// Sequences stored, not yet processed.
        /// </summary>
        private List<SAMAlignedSequence> sequences = null;

        /// <summary>
        /// Output queue for sequences which need to be written to a new filtered BAM file
        /// </summary>
        private Queue<List<SAMAlignedSequence>> bamOutputQueue = null;

        #endregion

        #region storage for sequence stats

        /// <summary>
        /// Average frequencies for each sequence, calculated per sample per cluster and averaged
        /// out to cluster
        /// </summary>
        private List<double[]> clustSeqFrequencies = new List<double[]>();

        /// <summary>
        /// The data that goes in a poisson graph, being <metric.CountAll, numClusters>
        /// </summary>
        Dictionary<int, int> graphDataAllReads = new Dictionary<int, int>();

        /// <summary>
        /// The data that goes in a poisson graph, being <metric.CountDistinct, numClusters>
        /// </summary>
        Dictionary<int, int> graphDataDistinctReads = new Dictionary<int, int>();

        /// <summary>
        /// The data that goes in a poisson graph, being <metric.CountSamples, numClusters>
        /// </summary>
        Dictionary<int, int> graphDataIndividualsCounts = new Dictionary<int, int>();

        /// <summary>
        /// ..
        /// average is rounded to the nearest whole number
        /// </summary>
        Dictionary<int, int> graphDataIndividualsTotalReadcounts = new Dictionary<int, int>();

        /// <summary>
        /// ..
        /// average is rounded to the nearest whole number
        /// </summary>
        Dictionary<int, int> graphDataIndividualsDistinctReadcounts = new Dictionary<int, int>();

        /// <summary>
        /// Percentage of the population that is represented in the cluster (i.e. for a population of 48 sampled, 
        /// 24 would be 0.5)
        /// </summary>
        //private List<double[]> clustPercentagePopulation = new List<double[]>();

        /// <summary>
        /// The (unnormalised) total number of reads in each cluster
        /// Can be used to determine if the cluster is above or below average size, once the size for every cluster is known
        /// </summary>
        private List<int> clustActualSize = new List<int>();

        #endregion

        #region IO

        /// <summary>
        /// Formatter used to write metric data to file (metrics for input file)
        /// </summary>
        private MetricFormatter formatterOriginalFile = null;

        /// <summary>
        /// Formatter used to write metric data to file (metrics for filtered 'good' sequences)
        /// </summary>
        private MetricFormatter formatterFilteredFile = null;

        /// <summary>
        /// Used to write read data back out to filtered BAM file
        /// </summary>
        private Stream bamStream = null;

        /// <summary>
        /// Used to write read data back out to filtered BAM file
        /// </summary>
        private BAMFormatter bamFormatter = null;

        #endregion

        #region fields from gui form
        /// <summary>
        /// Field from GUI form. Should a new filtered bam file be created
        /// </summary>
        private bool writeToFilteredBam = true;

        /// <summary>
        /// Field from GUI form. Should a metric file be created for the input file
        /// </summary>
        private bool writeClusterMetricOriginal = true;

        /// <summary>
        /// Field from GUI form. Should a metric file be created for the filtered output file
        /// </summary>
        private bool writeClusterMetricFiltered = true;

        /// <summary>
        /// Field from GUI form. Should an overview metric file be created for the input file
        /// </summary>
        private bool writeOverviewMetricOriginal = true;

        /// <summary>
        /// Field from GUI form. Should an overview metric file be created for the filtered output file
        /// </summary>
        private bool writeOverviewMetricFiltered = true;

        #endregion

        #region good/bad filter cutoffs

        /// <summary>
        /// Max allowed dirt
        /// </summary>
        private double dirtCutoff = 1;

        /// <summary>
        /// Max allowed 'g dirt'
        /// </summary>
        private double gDirtCutoff = double.MaxValue;

        /// <summary>
        /// Min allowed alignment qualtiy (as average per cluster)
        /// </summary>
        private double alignQualCutoff = 0;

        /// <summary>
        /// Min allowed read qualtiy (as average per cluster)
        /// </summary>
        private double readQualCutoff = 0;

        /// <summary>
        /// Min allowed population percentage within each cluster
        /// </summary>
        private double popPercent = 0;

        /// <summary>
        /// The number of samples to expect for the current input file
        /// </summary>
        private int numSamples = 0;

        #endregion

        #region running totals and averages

        /// <summary>
        /// Number of clusters parsed so far
        /// </summary>
        private int numClustersParsed = 0;

        /// <summary>
        /// Number of good clusters so far
        /// </summary>
        private int goodCount = 0;

        /// <summary>
        /// Max number of samples found in a cluster so far
        /// </summary>
        private int maxSampleCount = 0;

        /// <summary>
        /// Max mapping quality found in a cluster so far
        /// </summary>
        private double maxMapQ = 0;

        /// <summary>
        /// Max read quality found in a cluster so far
        /// </summary>
        private double maxReadQ = 0;

        /// <summary>
        /// Running total for dirt (all clusters)
        /// </summary>
        private double totalDirt = 0;

        /// <summary>
        /// Running total for mapping qualtiy (all clusters)
        /// </summary>
        private double totalMapQ = 0;

        /// <summary>
        /// Running total for read quality (all clusters)
        /// </summary>
        private double totalReadQ = 0;

        /// <summary>
        /// Running total for dirt (good clusters)
        /// </summary>
        private double totalDirtGood = 0;

        /// <summary>
        /// Running total for mapping qualtiy (good clusters)
        /// </summary>
        private double totalMapQGood = 0;

        /// <summary>
        /// Running total for read quality (good clusters)
        /// </summary>
        private double totalReadQGood = 0;

        /// <summary>
        /// Average dirt (so far) for all clusters
        /// </summary>
        private double averageDirt = 0;

        /// <summary>
        /// Average mapping quality (so far) for all clusters
        /// </summary>
        private double averageMapQ = 0;

        /// <summary>
        /// Average read qualtiy (so far) for all clusters
        /// </summary>
        private double averageReadQ = 0;

        /// <summary>
        /// Average dirt (so far) for good clusters
        /// </summary>
        private double averageDirtGood = 0;

        /// <summary>
        /// Average mapping quality (so far) for good clusters
        /// </summary>
        private double averageMapQGood = 0;

        /// <summary>
        /// Average read quality (so far) for good clusters
        /// </summary>
        private double averageReadQGood = 0;

        #endregion

        
        #endregion

        #region Constructors

        /// <summary>
        /// Default constructor.
        /// </summary>
        public ClusterMetricHandlerPloidulator() { }

        /// <summary>
        /// Non-default constructor, used to set the file name
        /// </summary>
        /// <param name="clusterFileName">Name of the file metric data is to be written to.</param>
        public ClusterMetricHandlerPloidulator(string clusterFileName, int ploidy, 
            double dirtCutoff, double alignQualCutoff, double readQualCutoff, double popPercent, int numSamples, bool? outputToFile)
            : this()
        {
            this.fileName = clusterFileName;
            this.expectedPloidy = ploidy;
            this.dirtCutoff = dirtCutoff;
            this.alignQualCutoff = alignQualCutoff;
            this.readQualCutoff = readQualCutoff;
            this.popPercent = popPercent;
            this.numSamples = numSamples;
            this.writeToFilteredBam = (outputToFile.HasValue) ? outputToFile.Value : false;
        }

        #endregion

        #region Properties

        /// <summary>
        /// Id of the cluster to which current sequences belong.
        /// </summary>
        public string CurrCluster { get { return currCluster; } }

        /// <summary>
        /// number of clusters which have been parsed
        /// </summary>
        public int ClusterCount { get { return clusterCount; } }

        /// <summary>
        /// The header for the bam input file
        /// </summary>
        public SAMAlignmentHeader InputHeader { get { return header; } set { header = value; } }

        #region filters

        /// <summary>
        /// Maximum allowed 'g dirt' cutoff
        /// </summary>
        public double GDirtCutoff { get { return gDirtCutoff; } set { gDirtCutoff = value; } }

        #endregion

        #region output files

        /// <summary>
        /// Field from GUI form. Should a new filtered bam file be created
        /// </summary>
        public bool WriteToFilteredBam { get { return writeToFilteredBam; } set { writeToFilteredBam = value; } }

        /// <summary>
        /// Field from GUI form. Should a metric file be created for the input file
        /// </summary>
        public bool WriteClusterMetricOriginal { get { return writeClusterMetricOriginal; } set { writeClusterMetricOriginal = value; } }
        
        /// <summary>
        /// Field from GUI form. Should a metric file be created for the filtered output file
        /// </summary>
        public bool WriteClusterMetricFiltered { get { return writeClusterMetricFiltered; } set { writeClusterMetricFiltered = value; } }
        
        /// <summary>
        /// Field from GUI form. Should an overview metric file be created for the input file
        /// </summary>
        public bool WriteOverviewMetricOriginal { get { return writeOverviewMetricOriginal; } set { writeOverviewMetricOriginal = value; } }
        
        /// <summary>
        /// Field from GUI form. Should an overview metric file be created for the filtered output file
        /// </summary>
        public bool WriteOverviewMetricFiltered { get { return writeOverviewMetricFiltered; } set { writeOverviewMetricFiltered = value; } }

        #endregion

        #region running totals and averages

        /// <summary>
        /// Number of good clusters so far
        /// </summary>
        public int GoodCount { get { return goodCount; } }

        /// <summary>
        /// Max number of samples found in a cluster so far
        /// </summary>
        public int MaxSampleCount { get { return maxSampleCount; } }

        /// <summary>
        /// Max mapping quality found in a cluster so far
        /// </summary>
        public double MaxMapQ { get { return maxMapQ; } }

        /// <summary>
        /// Max read quality found in a cluster so far
        /// </summary>
        public double MaxReadQ { get { return maxReadQ; } }

        /// <summary>
        /// Average dirt (so far) for all clusters
        /// </summary>
        public double AverageDirt { get { return averageDirt; } }

        /// <summary>
        /// Average mapping quality (so far) for all clusters
        /// </summary>
        public double AverageMapQ { get { return averageMapQ; } }

        /// <summary>
        /// Average read qualtiy (so far) for all clusters
        /// </summary>
        public double AverageReadQ { get { return averageReadQ; } }

        /// <summary>
        /// Average dirt (so far) for good clusters
        /// </summary>
        public double AverageDirtGood { get { return averageDirtGood; } }

        /// <summary>
        /// Average mapping quality (so far) for good clusters
        /// </summary>
        public double AverageMapQGood { get { return averageMapQGood; } }

        /// <summary>
        /// Average read quality (so far) for good clusters
        /// </summary>
        public double AverageReadQGood { get { return averageReadQGood; } }

        #endregion

        #endregion

        #region Delegates

        /// <summary>
        /// Delegate for performing cluster processing
        /// </summary>
        private delegate System.Delegate ClusterDelegate();

        #endregion

        #region Public Methods

        #region add

        /// <summary>
        /// Add a list of sequences. For each sequence, if the sequence belongs to the current cluster, 
        /// store it. If the sequence is part of a new cluster, process the current sequence cluster 
        /// then add the sequenceto a new cluster
        /// </summary>
        /// <param name="sequences">A list of sequences.</param>
        /// <returns>Always returns true</returns>
        public bool AddRange(IEnumerable<SAMAlignedSequence> sequences)
        {
            foreach (SAMAlignedSequence seq in sequences)
            {
                Add(seq); 
            }
            return true;
        }

        /// <summary>
        /// Add a sequence. If the sequence belongs to the current cluster, store it. If the sequence
        /// is part of a new cluster, process the current sequence cluster then add the sequence
        /// to a new cluster
        /// </summary>
        /// <param name="sequence">A sequence.</param>
        /// <returns>Returns true if the sequence could be added, false if the handler has been closed</returns>
        public bool Add(SAMAlignedSequence sequence)
        {
            if (sequences == null)
            {
                sequences = new List<SAMAlignedSequence>();
            }

            if (!finished)
            {
                string thisSeqCluster = sequence.RName; // Cluster the sequence we just added belongs to

                // This is the first sequence for the first cluster
                if (currCluster == null)
                {
                    currCluster = thisSeqCluster;
                }

                // This sequence belongs to a different cluster from the ones currently stored by this handler
                // (Process currently stored sequences before adding the new sequence)
                else if (!currCluster.Equals(thisSeqCluster)) 
                {
                    ++numClustersParsed; // mark off another cluster
                    ProcessSequences();
                    
                    currCluster = thisSeqCluster;
                    sequences = new List<SAMAlignedSequence>();
                }
                
                sequences.Add(sequence);
                return true;
            }

            // Processing of sequences should be finished but we are still outputting to the bam file
            // Or we are supposed to write to a bam file and the header and body files have not yet been merged
            // Wait for output to the bam file to complete
            else if ((!canWriteToBam && writeToFilteredBam) || (!bamFilesMerged && writeToFilteredBam))
            {
                while (!canWriteToBam || !bamFilesMerged)
                {
                    Thread.Sleep(20000); // sleep 20 seconds
                }
                return true;
            }

            // finished == true, bam file is writable and bam files have been merged (or no bam file was ever written to)
            // returning false indicates to calling process that no more sequences will be accepted
            else 
            {
                return false;
            }
        }

        #endregion

        #region modify behaviour

        /// <summary>
        /// Abort the handler when safe to do so (after the current cluster has finished processing)
        /// </summary>
        public void Abort()
        {
            aborted = true;
        }

        #endregion

        #region processing

        /// <summary>
        /// Calculate metric for a single cluser of sequences (all stored sequences), 
        /// and write metric data to file.
        /// </summary>
        public void ProcessSequences()
        {
            bool isGood = true;
            if (sequences != null && sequences.Count > 0 && !isComplete) // do the following only if there are sequences to be processed
            {
                clusterCount++;
                ClusterMetricPloidulator metric = new ClusterMetricPloidulator(expectedPloidy, numSamples);
                                
                // Initialise metric output file/s
                InitMetricOutputFiles();

                // Initialise bam output file/s
                InitBamOutputFiles();

                // Perform core metric calculations on cluster sequences
                metric.Calculate(sequences);
                isGood = GoodOrBad(metric);

                // Get haplotype information
                int numHap;
                if(isGood){
                    WritePhaseGenotypeInput(metric);
                    numHap = TestPhase(metric);
                    if(numHap > 4)
                    {
                        Console.WriteLine("Hap count means cluster is not good");
                        isGood = false;
                    }
                } else {
                    numHap = -1;
                }
                
                Console.Write("Cluster has "+numHap + " haplotypes");

                // Get statistics from the metric for this new cluster
                CreateSummaryArrays(metric);
                SetOverviewStats(metric, isGood);
                
                // Output sequences to metric file/s and/or filtered bam file
                WriteToMetricOutputFiles(metric, isGood);
                AddToOutputBamQueueOrDispose(metric, isGood);

                // If the bam file is not currently being written to, and there are sequences in the queue ready to be
                // written, launch a new thread to perform the writing to file
                if (writeToFilteredBam && canWriteToBam && bamOutputQueue.Count > 0)
                {
                    canWriteToBam = false;
                    ClusterDelegate runner = new ClusterDelegate(WriteToBam);
                    runner.BeginInvoke(null, null);
                }        
            }

            // Now that all processing has been performed for the current cluster, if the handler has been 
            // aborted, perform any final file outputs
            if(aborted)
            {
                SetComplete(false);
            }
        }

        #endregion

        #region completion
        /// <summary>
        /// Process any/all remaining sequences.
        /// </summary>
        public void SetComplete()
        {
            SetComplete(true);
        }

        /// <summary>
        /// Set complete but do not process any more sequences
        /// </summary>
        public void SetComplete(bool checkSequences)
        {
            if (!isComplete)
            {
                if (checkSequences)
                {
                    ProcessSequences();
                }

                PrintFinalCluserMetrics();

                if (writeToFilteredBam)
                {
                    concatFiles();
                }
                isComplete = true;
            }
            // Any other cleanup required that doesn't fit into Dispose()? todo
        }

        /// <summary>
        /// Disposes all formatters
        /// </summary>
        public void Dispose()
        {
            //SetComplete();
            // todo fixme i have probably done this wrong. there will be some other checking and stuff
            // that i also need to do or other stuff I need to clean up
            if (formatterOriginalFile != null)
            {
                formatterOriginalFile.Close();
            }
            if (formatterFilteredFile != null)
            {
                formatterFilteredFile.Close();
            }
            if (bamStream != null)
            {
                bamStream.Close();
            }
        }

#endregion

        #endregion

        #region Private Methods

        /// <summary>
        /// If metric is good, add it to the output queue (ready to be written to new BAM file)
        /// </summary>
        private void AddToOutputBamQueueOrDispose(ClusterMetricPloidulator metric, bool isGood)
        {
            if (writeToFilteredBam && isGood)
            {
                Console.WriteLine("Writing good cluster to file for # sequences: " + sequences.Count);

                // If the output queue has too many sequences in it, wait for the bam writer to catch up
                // and prevent memory fault
                if (bamOutputQueue.Count >= OUTPUT_QUEUE_SIZE)
                {
                    Console.WriteLine("Parser/processer thread must pause to let bam writer thread catch up");
                    while (bamOutputQueue.Count >= OUTPUT_QUEUE_SIZE / 2)
                    {
                        Thread.Sleep(10000); // sleep 10 seconds
                    }
                }

                // Add sequences to output file queue and output file header
                bamOutputQueue.Enqueue(sequences);
                AddToHeader(sequences[0]);
            }
            else
            {
                metric.Reset();
                metric = null;
                GC.Collect();
            }
        }

        /// <summary>
        /// Add a sequence to the filtered output file header
        /// </summary>
        private void AddToHeader(SAMAlignedSequence seq)
        {
            newHeader.ReferenceSequences.Add(new ReferenceSequenceInfo(seq.RName, 45)); // todo fixme this dataset only

            // for each good clust
            SAMRecordField sq = new SAMRecordField("SQ");
            sq.Tags.Add(new SAMRecordFieldTag("SN", seq.RName));
            sq.Tags.Add(new SAMRecordFieldTag("LN", "45")); // todo fixme this dataset only
            newHeader.RecordFields.Add(sq);
        }

        /// <summary>
        /// Write metrics to all per-cluster metric files
        /// </summary>
        private void WriteToMetricOutputFiles(ClusterMetricPloidulator metric, bool isGood)
        {
            if (writeClusterMetricOriginal)
            {
                formatterOriginalFile.Write(metric);
            }
            if (writeClusterMetricFiltered && isGood)
            {
                formatterFilteredFile.Write(metric);
            }
        }

        // should prob go within the metric instead
        private void WritePhaseGenotypeInput(ClusterMetricPloidulator metric)
        {
            int numIndividuals = metric.CountSamples;
            string locusType = metric.PhaseLoci; // S for a biallelic (SNP) locus; M for microsatellite, or other multi-allelic locus (eg tri-allelic SNP, or HLA allele).
            int numLoci = locusType.Length; // only the snps?
            string data = metric.PhaseData;
            //locusType = Regex.Split(locusType, "-").ToString();



            string lines = numIndividuals + "\r\n" +
                numLoci + "\r\n" +
                //"P 300 1313 1500 2023 5635\r\n"+ // position. this is optional and does not apply for sample dataset
                locusType + "\r\n" +
                data;

            //Console.Write(lines);

            System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "\\genotypes.inp");
            file.WriteLine(lines);
            file.Close();
        }



        private void WritePhaseGenotypeInputTemplate(ClusterMetricPloidulator metric)
        {
            int numIndividuals = 3;
            int numLoci = 5;
            string locusType = "MSSSM"; // S for a biallelic (SNP) locus; M for microsatellite, or other multi-allelic locus (eg tri-allelic SNP, or HLA allele).

            string indivA = "#1\r\n" +
                "12 1 0 1 3\r\n" +
                "11 0 1 0 3\r\n";

            string indivB = "#2\n12 1 1 1 2\r\n" +
                "12 0 0 0 3\r\n";

            string indivC = "#3\r\n" +
                "-1 ? 0 0 2\r\n" +
                "-1 ? 1 1 13";

            string lines = numIndividuals + "\r\n" +
                numLoci + "\r\n" +
                //"P 300 1313 1500 2023 5635\r\n"+ // position. this is optional and does not apply for sample dataset
                locusType + "\r\n" +
                indivA +
                indivB +
                indivC;

            System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "\\genotypes.inp");
            file.WriteLine(lines);
            file.Close();
        }

        /// <summary>
        /// Given a populated and calculated metric, determine based on handler's filter criteria whether
        /// that cluster is good or bad. Returns true for good, false for bad
        /// </summary>
        private bool GoodOrBad(ClusterMetricPloidulator tempMetric)
        {
            if (tempMetric.Dirt > dirtCutoff
                    || tempMetric.AlignmentQualities[0] < alignQualCutoff
                    || tempMetric.ReadQualities[0] < readQualCutoff
                    || tempMetric.PopulationPercentage < popPercent)
            {
                Console.WriteLine("BAD CLUSTER");
                tempMetric.Good = false;
            }
            else
            {
                Console.WriteLine("GOOD CLUSTER");
                ++goodCount;
                tempMetric.Good = true;
            }
            return tempMetric.Good;
        }

        /// <summary>
        /// Initialise all metric output files, if required and if not already initialised
        /// </summary>
        private void InitMetricOutputFiles()
        {
            if (writeClusterMetricOriginal && formatterOriginalFile == null)
            {
                formatterOriginalFile = new MetricFormatter(fileName + "\\orig.metr");
            }
            if (writeClusterMetricFiltered && formatterFilteredFile == null)
            {
                formatterFilteredFile = new MetricFormatter(fileName + "\\filtered.metr");
            }
        }

        /// <summary>
        /// Initialise bam output file, if required and if not already initialised
        /// </summary>
        private void InitBamOutputFiles()
        {
            if (writeToFilteredBam && bamStream == null && bamFormatter == null)
            {
                bamFormatter = new BAMFormatter();
                newHeader = new SAMAlignmentHeader();
                bamOutputQueue = new Queue<List<SAMAlignedSequence>>();

                // Create a new directory for output files if it does not already exist
                if (!Directory.Exists(fileName))
                {
                    Directory.CreateDirectory(fileName);
                }

                // Create the output file for filtered sequences
                string file = fileName + "\\sequences.bam";
                if(File.Exists(file))
                {
                    File.Delete(file);
                }
                bamStream = File.Create(file);
            }
        }

        /// <summary>
        /// Get various data arrays from metric, representing summary details for all clusters so far
        /// </summary>
        /// <param name="metric"></param>
        private void CreateSummaryArrays(ClusterMetricPloidulator metric) {
                    
            clustSeqFrequencies.Add(metric.ClustSeqFrequencies);
                    
            SetDictValueCounts(graphDataAllReads, metric.CountAll);
                    
            SetDictValueCounts(graphDataDistinctReads, metric.CountDistinct);
                    
            SetDictValueCounts(graphDataIndividualsCounts, metric.CountSamples);
                    
            int key = (int)Math.Round(metric.SampleReadCountsDistinct.Average(), 1);
            SetDictValueCounts(graphDataIndividualsDistinctReadcounts, key);
                    
            key = (int)Math.Round(metric.SampleReadCountsAll.Average(), 1);
            SetDictValueCounts(graphDataIndividualsTotalReadcounts, key);
        }

        /// <summary>
        /// Set or update overview stats
        /// (Totals are used to enable easy obtaining of average without iterating through lists to count each time)
        /// </summary>
        private void SetOverviewStats(ClusterMetricPloidulator metric, bool isGood)
        {
            // greatest number of samples found so far in any one cluster
            maxSampleCount = (metric.CountSamples > maxSampleCount) ? metric.CountSamples : maxSampleCount;

            // maximum quality value found so far in any one cluster
            maxMapQ = (metric.AlignmentQualities[0] > maxMapQ) ? metric.AlignmentQualities[0] : maxMapQ;
            maxReadQ = (metric.ReadQualities[0] > maxReadQ) ? metric.ReadQualities[0] : maxReadQ;

            // set totals for all clusters
            totalDirt += metric.Dirt;
            totalMapQ += metric.AlignmentQualities[0];
            totalReadQ += metric.ReadQualities[0];

            // set totals for good clusters
            totalDirtGood = (isGood) ? totalDirtGood + metric.Dirt : totalDirtGood;
            totalMapQGood += (isGood) ? metric.AlignmentQualities[0] : 0;
            totalReadQGood += (isGood) ? metric.ReadQualities[0] : 0;

            // running averages for all clusters
            averageDirt = Math.Round(totalDirt / (double)numClustersParsed, 2);
            averageMapQ = Math.Round(totalMapQ / (double)numClustersParsed, 2);
            averageReadQ = Math.Round(totalReadQ / (double)numClustersParsed, 2);

            // running averages for good clusters
            averageDirtGood = (goodCount > 0) ? Math.Round(totalDirtGood / (double)goodCount, 2) : 0;
            averageMapQGood = (goodCount > 0) ? Math.Round(totalMapQGood / (double)goodCount, 2) : 0;
            averageReadQGood = (goodCount > 0) ? Math.Round(totalReadQGood / (double)goodCount, 2) : 0;
        }

        /// <summary>
        /// returns the number of haplotypes
        /// </summary>
        private int TestPhase(ClusterMetricPloidulator metric)
        {
            Console.WriteLine("TOTAL G DIRT IS " + metric.TotalGDirt);
            if(metric.TotalGDirt == 0) // no alleles are outside the top 3
            {
                Console.WriteLine("Running PHASE to calculate haplotypes. Please wait...");
                ProcessStartInfo startInfo = new ProcessStartInfo();
                //startInfo.CreateNoWindow = false;
                startInfo.UseShellExecute = true;
                startInfo.FileName = "PHASE.exe";
                startInfo.ErrorDialog = false;
                startInfo.WindowStyle = ProcessWindowStyle.Hidden;
                startInfo.Arguments = "-d1 "+fileName+"\\genotypes.inp "+fileName+"\\haplotypes.out";

                try
                {
                    // Start the process and wait for it to finish
                    using (Process exeProcess = Process.Start(startInfo))
                    {
                        exeProcess.WaitForExit();
                    }
                }
                catch
                {
                    // Log error todo
                }
                Console.WriteLine("DONE");
                return GetPhaseNumHaplotypes(fileName + "\\haplotypes.out");
            }
            else
            {
                Console.WriteLine("Unable to get haplotypes, too many errors in data (gdirt " + metric.TotalGDirt + ")");
            }

            return -1;
        }

        private int GetPhaseNumHaplotypes(string file)
        {
            string hapIndex = "0";
            bool inHapBlock = false;
            try
            {
                using (StreamReader sr = new StreamReader(file))
                {
                    while (!sr.EndOfStream)
                    {
                        string thisStr = sr.ReadLine().Trim();
                        if (thisStr == "BEGIN LIST_SUMMARY")
                        {
                            inHapBlock = true;
                            thisStr = sr.ReadLine().Trim();
                        }
                        if(inHapBlock)
                        {
                            if (thisStr != "END LIST_SUMMARY")
                            {
                                hapIndex = thisStr.Split(' ')[0];
                            }
                            else
                            {
                                Console.WriteLine("Best haplotype count for this cluster is " + hapIndex);
                                return Convert.ToInt32(hapIndex);
                            }
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("The file could not be read:");
                Console.WriteLine(e.Message);
            }
            return -1;
        }

        /// <summary>
        /// Concatenate together the header and sequence bam output files
        /// </summary>
        private void concatFiles()
        {
            Console.WriteLine("Concatenating bam header and sequence output files");
            ProcessStartInfo startInfo = new ProcessStartInfo();

            startInfo.CreateNoWindow = true;
            startInfo.FileName = "cmd.exe";
            startInfo.Arguments = "/C copy /b " + fileName + "\\header.bam+" + fileName + "\\sequences.bam " + fileName + "\\filtered.bam /y";
            Console.WriteLine("Executing cmd "+startInfo.Arguments);
            
            startInfo.UseShellExecute = false;
            startInfo.RedirectStandardOutput = true;
            

            using (Process proc = Process.Start(startInfo))
            {
                using (StreamReader reader = proc.StandardOutput)
                {
                    string result = reader.ReadToEnd();
                    Console.WriteLine(result);
                }
                proc.WaitForExit();

            }

            bamFilesMerged = true;
            Console.WriteLine("Finished concatenating");

            //Process.Start(OUTPUT_FOLDER);
        }
        

        /// <summary>
        /// For the given key-count dictionary, if it has key, increment the value to reflect that it contains
        /// another instance. If it does not yet have that key, add a first instance with that key reference
        /// Dict values are rounded to the nearest whole number (int)
        /// </summary>
        private void SetDictValueCounts(Dictionary<int, int> dict, int key)
        {
            if (dict.ContainsKey(key))
            {
                ++dict[key];
            }
            else
            {
                dict[key] = 1;
            }
        }
        
        /*
         I haven't cleaned this method up because it needs to be significantly rewritten
         */
        private void PrintFinalCluserMetrics()
        {
            if (!finished)
            {

                // writes to a separate file for each significant all-clusters metric

                /*Gives us the actual, unnormalised size of each cluster (independent of which individuals the reads come from)
                 * Can be used to determine if each is above/below average size and/or the degree of distribution/variation
                 * 
                 * <numReadsAll, numClusters>
                 * <numReadsDistinct, numClusters>
                 * 
                 * 
                 * 1. normalise each
                 */
                #region clusters_avg_size
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_allreadcounts.tsv"))
                {
                    file.WriteLine("num_reads\tnum_clusters");
                    foreach (KeyValuePair<int, int> dat in graphDataAllReads)
                    {
                        file.WriteLine(dat.Key + "\t" + dat.Value);
                    }
                }

                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_distinctreadcounts.tsv"))
                {
                    file.WriteLine("num_reads\tnum_clusters");
                    foreach (KeyValuePair<int, int> dat in graphDataDistinctReads)
                    {
                        file.WriteLine(dat.Key + "\t" + dat.Value);
                    }
                }
                #endregion


                /*
             * How many individuals are represented in each cluster? Also as Poisson-graphable type data
             * e.g. proportion of clusters that have all population represented Vs proportion that have half population represnted
             * <numIndividuals, numClusters>
             */
                #region cluster_individual_representation

                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_numindividualscounts.tsv"))
                {
                    file.WriteLine("num_reads\tnum_clusters");
                    foreach (KeyValuePair<int, int> dat in graphDataIndividualsCounts)
                    {
                        file.WriteLine(dat.Key + "\t" + dat.Value);
                    }
                }

                // even if every individual is represented in the cluster, is there a vast difference between how they are represented?
                // something to work out how many reads on average each individual has (read count per indiv averaged across cluster instead of just
                // read count per cluster)
                // average total readcount by individual
                //graphDataIndividualsTotalReadcounts
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_byindividualsrtotaleadcounts.tsv"))
                {
                    file.WriteLine("num_reads\tnum_clusters");
                    foreach (KeyValuePair<int, int> dat in graphDataIndividualsTotalReadcounts)
                    {
                        file.WriteLine(dat.Key + "\t" + dat.Value);
                    }
                }

                // average distinct readcount by individual
                //graphDataIndividualsDistinctReadcounts
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_byindividualsrdistincteadcounts.tsv"))
                {
                    file.WriteLine("num_reads\tnum_clusters");
                    foreach (KeyValuePair<int, int> dat in graphDataIndividualsDistinctReadcounts)
                    {
                        file.WriteLine(dat.Key + "\t" + dat.Value);
                    }
                }

                #endregion

                /*
             * The frequency distributions of the reads in each cluster on a per-individual basis
             * Normalised as percentages from smaller to larger (e.g. piechart-type data)
             * 
             * 
             */
                #region frequency_distributions
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_frequencies.tsv"))
                {
                    file.WriteLine("frequencies");
                    foreach (double[] dat in clustSeqFrequencies)
                    {
                        file.WriteLine(String.Join("\t", dat));
                    }
                }
                #endregion

                Console.WriteLine("finished");

                finished = true;
            }
            while (!canWriteToBam) // background thread is currently writing
            {
                Thread.Sleep(5000); // sleep 5 seconds
            }
            if(writeToFilteredBam){
                canWriteToBam = false;
                string f = fileName + "\\header.bam";
                if (File.Exists(f))
                {
                    File.Delete(f);
                }
                using (Stream writer = File.OpenWrite(f))
                {
                    new BAMFormatter().WriteHeader(newHeader, writer);
                }
                //Thread.Sleep(5000); // sleep 5 seconds
                canWriteToBam = true;
            }
        }

        /// <summary>
        /// Write sequences from output queue into filtered bam output file
        /// </summary>
        /// <returns></returns>
        private System.Delegate WriteToBam()
        {
            while (bamOutputQueue.Count > 0)
            {
                List<SAMAlignedSequence> seqs = bamOutputQueue.Dequeue();
                foreach (SAMAlignedSequence seq in seqs)
                {
                    bamFormatter.WriteAlignedSequence(header, seq, bamStream);
                }
                seqs.Clear();
                seqs = null;
                GC.Collect();
                
            }
            // signal to next thread runner that it can now process sequences from the queue
            canWriteToBam = true;
            return null;
        }

        #endregion


    }
}
