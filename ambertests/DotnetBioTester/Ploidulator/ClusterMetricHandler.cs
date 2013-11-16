using Bio.IO.BAM;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading;

namespace Ploidulator
{
    /// <summary>
    /// An ClusterMetricHandler receives a set of SAMAlignedSequences (as a list or one by one)
    /// and passes them as "clusters" to a ClusterMetric, which performs a number of calculations
    /// on them to calculate the accuracy of the cluster.
    /// 
    /// All determination as to whether a cluster is 'good' or 'bad' is managed by ClusterMetricHandlerPloidulator
    /// </summary>
    public class ClusterMetricHandler : IMetricHandler
    {
        private static CultureInfo ci = new CultureInfo("en-AU");

        // for debugging purposes haplotyping can be disabled (faster calculation of all other metrics)
        private bool haplotypingEnabled = true;

        #region Private Static Fields
        private static int OUTPUT_QUEUE_SIZE = 7;                       // max number of clusters that can be stored in the
                                                                        // output queue (to prevent too many being held in memory)
        private static string OUTPUT_DIRECTORY = "orig";                // output directory for metrics relating to original file
        private static string FILTERED_OUTPUT_DIRECTORY = "filtered";   // output directory for metrics relating to filtered file
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
        /// Id of the cluster to which current sequences belong
        /// </summary>
        private string currentClusterId;

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
        /// Whether the handler has been aborted
        /// </summary>
        private bool aborted = false;

        /// <summary>
        /// Whether the handler has done its completion cleanup 
        /// </summary>
        private bool isComplete = false;

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

        /// <summary>
        /// This value is true if haplotypes.txt should be constructed
        /// </summary>
        private bool writeHaplotypesFile = false;

        /// <summary>
        /// This value is true if genotypes.txt should be constructed
        /// </summary>
        private bool writeGenotypesFile = false;

        #endregion

        #region storage for sequences

        /// <summary>
        /// Sequences stored for the current cluster, not yet processed.
        /// </summary>
        private Collection<SAMAlignedSequence> allSequences = null;

        /// <summary>
        /// Output queue for sequences which need to be written to a new filtered BAM file
        /// </summary>
        private Queue<Collection<SAMAlignedSequence>> bamOutputQueue = null;

        #endregion

        #region storage for sequence stats

        /// <summary>
        /// Average frequencies for each sequence, calculated per sample per cluster and averaged
        /// out to cluster
        /// </summary>
        private List<Collection<double>> clustSeqFrequencies = new List<Collection<double>>();

        /// <summary>
        /// A single frequency distribution, in the same format as one row from clustSeqFrequencies, intended to give
        /// an overview of frequency distribution across all clusters
        /// </summary>
        private Collection<double> clusterSequenceFrequenciesOverview = new Collection<double>();

        /// <summary>
        /// A single frequency distribution, in the same format as one row from clustSeqFrequencies, intended to give
        /// an overview of frequency distribution across all good clusters
        /// </summary>
        private Collection<double> clusterSequenceFrequenciesOverviewGood = new Collection<double>();

        /// <summary>
        /// Average frequencies for each good sequence, calculated per sample per cluster and averaged
        /// out to cluster
        /// </summary>
        private List<Collection<double>> clustSeqFrequenciesGood = new List<Collection<double>>();

        /// <summary>
        /// The data that goes in a poisson graph, being <metric.CountAll, numClusters>
        /// </summary>
        Dictionary<int, int> graphDataAllReads = new Dictionary<int, int>();

        /// <summary>
        /// The data that goes in a poisson graph, being <metric.CountAll, numClusters>
        /// </summary>
        Dictionary<int, int> graphDataAllReadsGood = new Dictionary<int, int>();

        /// <summary>
        /// The data that goes in a poisson graph, being <metric.CountDistinct, numClusters>
        /// </summary>
        Dictionary<int, int> graphDataDistinctReads = new Dictionary<int, int>();

        /// <summary>
        /// The data that goes in a poisson graph, being <metric.CountDistinct, numClusters>
        /// </summary>
        Dictionary<int, int> graphDataDistinctReadsGood = new Dictionary<int, int>();

        /// <summary>
        /// The data that goes in a poisson graph, being <metric.CountSamples, numClusters>
        /// </summary>
        Dictionary<int, int> graphDataIndividualsCounts = new Dictionary<int, int>();

        /// <summary>
        /// The data that goes in a poisson graph, being <metric.CountSamples, numClusters>
        /// </summary>
        Dictionary<int, int> graphDataIndividualsCountsGood = new Dictionary<int, int>();

        /// <summary>
        /// ..
        /// average is rounded to the nearest whole number
        /// </summary>
        Dictionary<int, int> graphDataIndividualsTotalReadcounts = new Dictionary<int, int>();

        /// <summary>
        /// ..
        /// average is rounded to the nearest whole number
        /// </summary>
        Dictionary<int, int> graphDataIndividualsTotalReadcountsGood = new Dictionary<int, int>();

        /// <summary>
        /// ..
        /// average is rounded to the nearest whole number
        /// </summary>
        Dictionary<int, int> graphDataIndividualsDistinctReadcounts = new Dictionary<int, int>();

        /// <summary>
        /// ..
        /// average is rounded to the nearest whole number
        /// </summary>
        Dictionary<int, int> graphDataIndividualsDistinctReadcountsGood = new Dictionary<int, int>();


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
        /// Used to write data back out to filtered BAM file
        /// </summary>
        private Stream bamStream = null;

        /// <summary>
        /// Used to write data back out to filtered BAM file
        /// </summary>
        private BAMFormatter bamFormatter = null;

        /// <summary>
        /// Used to write data to genotypes.txt
        /// </summary>
        private StreamWriter genotypesStream = null;

        /// <summary>
        /// Used to write data to haplotypes.txt
        /// </summary>
        private StreamWriter haplotypesStream = null;

        #endregion

        #region fields from gui form

        /// <summary>
        /// Field from GUI form. Should a new filtered bam file be created
        /// </summary>
        private bool writeToFilteredBam = true;

        /// <summary>
        /// Field from GUI form. Should every cluster be haplotyped or only ones which are probably good. False by default
        /// </summary>
        private bool onlyHaplotypeGood = false;

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
        /// Max allowed dirt. Default 1.
        /// </summary>
        private double dirtCutoff = 1;

        /// <summary>
        /// Max allowed haplotypes in a cluster
        /// </summary>
        private int hapMaxCutoff = int.MaxValue;

        /// <summary>
        /// Max allowed ploidy disagreement in a cluster
        /// </summary>
        private double ploidyDisagreementCutoff = int.MaxValue;

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
        private double populationPercentageCutoff = 0;

        /// <summary>
        /// The number of samples to expect for the current input file
        /// </summary>
        private int numSamples = 0;

        #endregion

        #region running totals and averages

        /// <summary>
        /// Number of clusters parsed so far
        /// </summary>
        private int numberClustersParsed = 0;

        /// <summary>
        /// Number of good clusters so far
        /// </summary>
        private int goodCount = 0;

        /// <summary>
        /// Max number of samples found in a cluster so far
        /// </summary>
        private int maxSampleCount = 0;

        /// <summary>
        /// Total number of reads
        /// </summary>
        private int readCountTotal = 0;

        /// <summary>
        /// Total number of 'good' reads
        /// </summary>
        private int readCountGood = 0;

        /// <summary>
        /// Number of distinct reads
        /// </summary>
        private int readCountDistinctTotal = 0;

        /// <summary>
        /// Number of 'good' distinct reads
        /// </summary>
        private int readCountDistinctGood = 0;

        /// <summary>
        /// Max mapping quality found in a cluster so far
        /// </summary>
        private double maxMapQuality = 0;

        /// <summary>
        /// Max read quality found in a cluster so far
        /// </summary>
        private double maxReadQuality = 0;

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
        /// Constructor, used to set the file name, and the ploidy level and number of individuals in that file
        /// </summary>
        /// <param name="clusterFileName">Name of the input file.</param>
        /// <param name="ploidy">Expected ploidy level of the organism in the input file.</param>
        /// <param name="numberSamples">Number of samples present in the input file.</param>
        public ClusterMetricHandler(string clusterFileName, int ploidy, int numberSamples)
        {
            this.fileName = clusterFileName;
            this.expectedPloidy = ploidy;
            this.numSamples = numberSamples;

            // Create a new output directory with name of input file, if it does not already exist
            if (!Directory.Exists(fileName))
            {
                Directory.CreateDirectory(fileName);
            }
        }

        #endregion

        #region Properties

        /// <summary>
        /// Get the name of the input file/directory (excluding file extension).
        /// </summary>
        public string FileName { get { return fileName; } }

        /// <summary>
        /// Get the id of the cluster to which current sequences belong.
        /// </summary>
        public string CurrentClusterId { get { return currentClusterId; } }

        /// <summary>
        /// Get the number of clusters which have been parsed (so far).
        /// </summary>
        public int ClusterCount { get { return clusterCount; } }

        /// <summary>
        /// Get the header for the bam input file.
        /// </summary>
        public SAMAlignmentHeader InputHeader { get { return header; } set { header = value; } }

        #region filters

        /// <summary>
        /// Get or set the maximum number haplotypes allowed in a good cluster
        /// </summary>
        public int HaplotypesMaxCutoff { get { return hapMaxCutoff; } set { hapMaxCutoff = value; } }

        /// <summary>
        /// Get or set the maximum allowed ploidy disagreement in a good cluster
        /// </summary>
        public double PloidyDisagreementCutoff { get { return ploidyDisagreementCutoff; } set { ploidyDisagreementCutoff = value; } }

        /// <summary>
        /// Get or set the maximum dirt allowed in a good cluster. Default 1.
        /// </summary>
        public double DirtCutoff { get { return dirtCutoff; } set { dirtCutoff = value; } }

        /// <summary>
        /// Get or set the minimum alignment qualtiy allowed in a good cluster (as average per cluster).
        /// </summary>
        public double AlignmentQualityCutoff { get { return alignQualCutoff; } set { alignQualCutoff = value; } }

        /// <summary>
        /// Get or set the minimum read quality allowed in a good cluster (as average per cluster).
        /// </summary>
        public double ReadQualityCutoff { get { return readQualCutoff; } set { readQualCutoff = value; } }

        /// <summary>
        /// Get or set the minimum population percentage allowed within each good cluster.
        /// </summary>
        public double PopulationPercentageCutoff { get { return populationPercentageCutoff; } set { populationPercentageCutoff = value; } }

        #endregion

        #region output files

        /// <summary>
        /// Get or set a flag to indicate whether a new filtered bam file be created.
        /// </summary>
        public bool WriteToFilteredBam { get { return writeToFilteredBam; } set { writeToFilteredBam = value; } }

        /// <summary>
        /// Get or set a flag to indicate whether haplotypes should be calculated for all clusters or only those which are probably good based on other criteria
        /// </summary>
        public bool OnlyHaplotypeGood { get { return onlyHaplotypeGood; } set { onlyHaplotypeGood = value; } }

        /// <summary>
        /// Get or set a flag to indicate whether a metric file should be created for the input file
        /// </summary>
        public bool WriteClusterMetricOriginal { get { return writeClusterMetricOriginal; } set { writeClusterMetricOriginal = value; } }
        
        /// <summary>
        /// Get or set a flag to indicate whether a metric file should be created for the filtered output file
        /// </summary>
        public bool WriteClusterMetricFiltered { get { return writeClusterMetricFiltered; } set { writeClusterMetricFiltered = value; } }
        
        /// <summary>
        /// Get or set a flag to indicate whether an overview metric file should be created for the input file
        /// </summary>
        public bool WriteOverviewMetricOriginal { get { return writeOverviewMetricOriginal; } set { writeOverviewMetricOriginal = value; } }
        
        /// <summary>
        /// Get or set a flag to indicate whether an overview metric file should be created for the filtered output file
        /// </summary>
        public bool WriteOverviewMetricFiltered { get { return writeOverviewMetricFiltered; } set { writeOverviewMetricFiltered = value; } }

        /// <summary>
        /// Get or set a flag to indicate whether genotypes.txt should be constructed
        /// </summary>
        public bool WriteGenotypesFile { get { return writeGenotypesFile; } set { writeGenotypesFile = value; } }

        /// <summary>
        /// Get or set a flag to indicate whether haplotypes.txt should be constructed
        /// </summary>
        public bool WriteHaplotypesFile { get { return writeHaplotypesFile; } set { writeHaplotypesFile = value; } }

        #endregion

        #region running totals and averages

        /// <summary>
        /// Get the number of good clusters so far.
        /// </summary>
        public int GoodCount { get { return goodCount; } }

        /// <summary>
        /// Get the total number of reads (so far).
        /// </summary>
        public int ReadCountTotal { get { return readCountTotal; } }

        /// <summary>
        /// Get the total number of 'good' reads (so far).
        /// </summary>
        public int ReadCountGood { get { return readCountGood; } }

        /// <summary>
        /// Get the total number of distinct reads (so far).
        /// </summary>
        public int ReadCountDistinctTotal { get { return readCountDistinctTotal; } }
        
        /// <summary>
        /// Get the total number of 'good' distinct reads (so far).
        /// </summary>
        public int ReadCountDistinctGood { get { return readCountDistinctGood; } }
        
        /// <summary>
        /// Get the number of clusters parsed so far.
        /// </summary>
        public int NumberClustersParsed { get { return numberClustersParsed; } }

        /// <summary>
        /// Get the maximum number of samples found in a cluster so far.
        /// </summary>
        public int MaxSampleCount { get { return maxSampleCount; } }

        /// <summary>
        /// Get the maximum alignment quality found in a cluster so far.
        /// </summary>
        public double MaxAlignmentQuality { get { return maxMapQuality; } }

        /// <summary>
        /// Get the maximum read quality found in a cluster so far.
        /// </summary>
        public double MaxReadQuality { get { return maxReadQuality; } }

        /// <summary>
        /// Get the average dirt (so far) for all clusters.
        /// </summary>
        public double AverageDirt { get { return averageDirt; } }

        /// <summary>
        /// Get the average alignment quality (so far) for all clusters.
        /// </summary>
        public double AverageMapQ { get { return averageMapQ; } }

        /// <summary>
        /// Get the average read qualtiy (so far) for all clusters.
        /// </summary>
        public double AverageReadQ { get { return averageReadQ; } }

        /// <summary>
        /// Get the average dirt (so far) for GOOD clusters.
        /// </summary>
        public double AverageDirtGood { get { return averageDirtGood; } }

        /// <summary>
        /// Get the average mapping quality (so far) for GOOD clusters.
        /// </summary>
        public double AverageMapQGood { get { return averageMapQGood; } }

        /// <summary>
        /// Get the average read quality (so far) for GOOD clusters.
        /// </summary>
        public double AverageReadQGood { get { return averageReadQGood; } }

        /// <summary>
        /// Get an overview of sequence frequency distribution across all clusters.
        /// </summary>
        public Collection<double> ClusterSequenceFrequenciesOverview { get { return clusterSequenceFrequenciesOverview; } }

        /// <summary>
        /// Get an overview of sequence frequency distribution across all GOOD clusters.
        /// </summary>
        public Collection<double> ClusterSequenceFrequenciesOverviewGood { get { return clusterSequenceFrequenciesOverviewGood; } }

        /// <summary>
        /// Get data for the population distribution graph (all clusters).
        /// </summary>
        public Dictionary<int,int> GraphDataIndividualsCounts 
        { 
            get 
            {
                return (from datum in graphDataIndividualsCounts orderby datum.Key descending select datum)
                    .ToDictionary(pair => (int)(pair.Key / (double)numSamples * 100), pair => (int)((double)pair.Value 
                        / (double)numberClustersParsed * (double)100));
            } 
        }

        /// <summary>
        /// Get data for the population distribution graph (good clusters).
        /// </summary>
        public Dictionary<int, int> GraphDataIndividualsCountsGood
        {
            get
            {
                return (from datum in graphDataIndividualsCountsGood orderby datum.Key descending select datum)
                    .ToDictionary(pair => (int)(pair.Key / (double)numSamples * 100), pair => (int)(pair.Value / (double)goodCount * 100));
            }
        
        }
        
        /// <summary>
        /// Get data for the distinct read distribution graph (all clusters).
        /// </summary>
        public Dictionary<int, int> GraphDataDistinctReads
        {
            get
            {
                return (from datum in graphDataDistinctReads orderby datum.Value descending select datum)
                        .ToDictionary(pair => pair.Key, pair => pair.Value);
            }
        }

        /// <summary>
        /// Get data for the distinct read distribution graph (GOOD clusters).
        /// </summary>
        public Dictionary<int, int> GraphDataDistinctReadsGood
        {
            get
            {
                return (from datum in graphDataDistinctReadsGood orderby datum.Value descending select datum)
                        .ToDictionary(pair => pair.Key, pair => pair.Value);
            }
        }

        /// <summary>
        /// Get the actual minimum distinct read count (all clusters), or Int32.MinValue if read counts are not known.
        /// </summary>
        public int MinDistinctReadCount
        {
            get { return (graphDataDistinctReads.Count > 0) ? graphDataDistinctReads.Keys.Min() : Int32.MinValue; }
        }

        /// <summary>
        /// Get the actual maximum distinct read count (all clusters), or Int32.MaxValue if read counts are not known.
        /// </summary>
        public int MaxDistinctReadCount
        {
            get { return (graphDataDistinctReads.Count > 0) ? graphDataDistinctReads.Keys.Max() : Int32.MaxValue; }
        }

        #endregion

        #endregion

        #region Delegates (background thread runners)

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
        /// <returns>Always returns true.</returns>
        public bool AddRange(IEnumerable<SAMAlignedSequence> sequences)
        {
            if (sequences == null || sequences.Count() == 0)
            {
                return true;
            }
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
        /// <returns>Returns true if the sequence could be added, false if the handler has been closed.</returns>
        public bool Add(SAMAlignedSequence sequence)
        {
            if (allSequences == null)
            {
                allSequences = new Collection<SAMAlignedSequence>();
            }
            if(sequence == null)
            {
                return true;
            }

            if (!finished)
            {
                string thisSeqCluster = sequence.RName; // Cluster the sequence we just added belongs to

                // This is the first sequence for the first cluster
                if (currentClusterId == null)
                {
                    currentClusterId = thisSeqCluster;
                }

                // This sequence belongs to a different cluster from the ones currently stored by this handler
                // (Process currently stored sequences before adding the new sequence)
                else if (!currentClusterId.Equals(thisSeqCluster)) 
                {
                    ++numberClustersParsed; // mark off another cluster
                    ProcessSequences();
                    
                    currentClusterId = thisSeqCluster;
                    allSequences = new Collection<SAMAlignedSequence>();
                }
                
                allSequences.Add(sequence);
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

        #region modify behaviour or end process

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
        /// and write metric data to file/s if required.
        /// </summary>
        public void ProcessSequences()
        {
            bool isGood = true;
            if (allSequences != null && allSequences.Count > 0 && !isComplete) // do the following only if there are sequences to be processed
            {
                clusterCount++;
                ClusterMetric metric = new ClusterMetric(expectedPloidy, numSamples);
                                
                // Initialise metric output file/s
                InitMetricOutputFiles();

                // Initialise bam output file/s
                InitBamOutputFiles();

                // Perform core metric calculations on cluster sequences
                metric.Calculate(allSequences);
                isGood = GoodOrBad(metric);

                // Get haplotype information
                if(haplotypingEnabled && expectedPloidy == 2)
                {
                    GetHaplotypeInfo(ref metric, ref isGood);
                }
                if(isGood) { ++goodCount; }
                Console.WriteLine(metric.ToString() + "\t" + (isGood ? Properties.Resources.GOOD_CLUSTER : Properties.Resources.BAD_CLUSTER));

                // Get statistics from the metric for this new cluster
                CreateSummaryArrays(metric, isGood);
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
        /// Set complete and process all remaining sequences.
        /// </summary>
        public void SetComplete()
        {
            SetComplete(true);
        }

        /// <summary>
        /// Set complete.
        /// </summary>
        /// <param name="checkSequences">Boolean value indicates whether stored sequences should be processed
        /// before setting complete.</param>
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
                DeletePhaseFiles();
                if (goodCount < 1)
                {
                    DeleteFilteredMetricFile();
                }
                isComplete = true;
            }
        }

        /// <summary>
        /// Dispose resources.
        /// </summary>
        /// <param name="includeManagedResources">Boolean value indicates whether managed resources should also be disposed.</param>
        protected virtual void Dispose(bool includeManagedResources)
        {
            if(includeManagedResources){
                if (formatterOriginalFile != null)
                {
                    formatterOriginalFile.Close();
                    formatterOriginalFile = null;
                }
                if (formatterFilteredFile != null)
                {
                    formatterFilteredFile.Close();
                    formatterFilteredFile = null;
                }
                if (bamStream != null)
                {
                    bamStream.Close();
                    bamStream = null;
                }
                if (genotypesStream != null)
                {
                    genotypesStream.Close();
                    genotypesStream = null;
                }
                if (haplotypesStream != null)
                {
                    haplotypesStream.Close();
                    haplotypesStream = null;
                }
            }
        }

        /// <summary>
        /// Disposes all formatters.
        /// </summary>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        #endregion

        #endregion

        #region Private Methods

        #region queue sequences for output

        /// <summary>
        /// If metric is good, add it to the output queue (ready to be written to new BAM file)
        /// </summary>
        private void AddToOutputBamQueueOrDispose(ClusterMetric metric, bool isGood)
        {
            if (writeToFilteredBam && isGood)
            {
                // If the output queue has too many sequences in it, wait for the bam writer to catch up
                // and prevent memory fault
                if (bamOutputQueue.Count >= OUTPUT_QUEUE_SIZE)
                {
                    Console.WriteLine(Properties.Resources.PROCESSING_THREAD_PAUSED);
                    while (bamOutputQueue.Count >= OUTPUT_QUEUE_SIZE / 2)
                    {
                        Thread.Sleep(10000); // sleep 10 seconds
                    }
                }

                // Add sequences to output file queue and output file header
                bamOutputQueue.Enqueue(allSequences);
                AddToHeader(allSequences[0]);
            }
            else
            {
                metric.Reset();
                metric = null;
            }
        }

        #endregion

        #region phase

        /// <summary>
        /// If a cluster is good, runs phase to get haplotype count for that cluster. If the number of haplotypes
        /// found is > hapMaxCutoff, sets isGood to false. The number of haplotypes found is also passed
        /// to the metric
        /// </summary>
        private void GetHaplotypeInfo(ref ClusterMetric metric, ref bool isGood)
        {
            int numHap;
            if (!onlyHaplotypeGood || onlyHaplotypeGood && isGood)
            {
                WritePhaseGenotypeInput(metric);
                numHap = CalculateClusterHaplotypes();
                if (numHap > hapMaxCutoff)
                {
                    Console.WriteLine(Properties.Resources.HAPLOTYPE_COUNT_BAD + numHap);
                    isGood = false;
                    metric.Good = false;
                }
            }
            else
            {
                numHap = -1;
            }
            metric.NumberOfHaplotypes = numHap;
        }

        /// <summary>
        /// Returns the optimum number of haplotypes detected in the cluster, or -1 if G dirt was too high so
        /// process was not run. Phase output files will also be saved to the output directory
        /// </summary>
        private int CalculateClusterHaplotypes()
        {
            ProcessStartInfo startInfo = new ProcessStartInfo();
            startInfo.UseShellExecute = true; // set to false to display output in console
            startInfo.FileName = "PHASE.exe";
            startInfo.ErrorDialog = false;
            startInfo.WindowStyle = ProcessWindowStyle.Hidden;
            startInfo.Arguments = "-d1 " + fileName + "\\genotypes.inp " + fileName + "\\haplotypes.out";

            try
            {
                using (Process exeProcess = Process.Start(startInfo))
                {
                    exeProcess.WaitForExit();
                }
            }
            catch(Exception e)
            {
                Console.WriteLine(Properties.Resources.PHASE_ERROR_EXECUTING);
                Console.WriteLine(e.Message);
            }
            if(writeHaplotypesFile)
            {
                AddHaplotypesToMasterFile();
            }

            return GetPhaseNumHaplotypes(fileName + "\\haplotypes.out");
        }

        /// <summary>
        /// Copy the contents of haplotypes.out to the end of haplotypes.txt
        /// TODO this should really be coupled in with GetPhaseNumHaplotypes to prevent reading the file twice
        /// </summary>
        private void AddHaplotypesToMasterFile()
        {
            if (haplotypesStream == null)
            {
                haplotypesStream = new System.IO.StreamWriter(fileName + "\\haplotypes.txt");
            }
            using (StreamReader sr = new StreamReader(fileName + "\\haplotypes.out"))
            {
                while (!sr.EndOfStream)
                {
                    haplotypesStream.WriteLine(sr.ReadLine());
                }
            }
        }

        /// <summary>
        /// Parse the primary PHASE output file to extract the number of haplotypes. Returns the optimum number of 
        /// haplotypes detected in the cluster, or -1 if the file could not be read
        /// </summary>
        private static int GetPhaseNumHaplotypes(string file)
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
                        if (inHapBlock)
                        {
                            if (thisStr != "END LIST_SUMMARY")
                            {
                                hapIndex = thisStr.Split(' ')[0];
                            }
                            else
                            {
                                return Convert.ToInt32(hapIndex, ci);
                            }
                        }
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(Properties.Resources.FILE_NOT_READABLE);
                Console.WriteLine(e.Message);
            }
            return -1;
        }

        /// <summary>
        /// Write PHASE input data for metric to file
        /// </summary>
        /// <param name="metric"></param>
        private void WritePhaseGenotypeInput(ClusterMetric metric)
        {
            int numIndividuals = metric.CountSamples;
            string locusType = metric.PhaseLoci; // S for a biallelic (SNP) locus; M for microsatellite, or other multi-allelic locus (eg tri-allelic SNP, or HLA allele).
            int numLoci = locusType.Length; 
            string data = metric.PhaseData;

            string lines = numIndividuals + "\r\n" +
                numLoci + "\r\n" +
                //"P 300 1313 1500 2023 5635\r\n"+ // position. this is optional and does not apply for sample dataset
                locusType + "\r\n" +
                data;

            // While the files cannot be deleted, either the previous phase.exe process is still writing to/reading from
            // one of these files, or the user has one open
            while (!DeletePhaseFiles())
            {
                Console.WriteLine(Properties.Resources.PHASE_ERROR_DELETE);
                Thread.Sleep(20000); // sleep 20 seconds before trying again
            }
            
            
            using (StreamWriter file = new System.IO.StreamWriter(fileName + "\\genotypes.inp"))
            {
                file.WriteLine(lines);
                if(writeGenotypesFile)
                {
                    if (genotypesStream == null)
                    {
                        genotypesStream = new System.IO.StreamWriter(fileName + "\\genotypes.txt");
                    }
                    genotypesStream.WriteLine(lines);
                }
            }
        }

        /// <summary>
        /// Deletes all PHASE output AND INPUT files.
        /// </summary>
        private bool DeletePhaseFiles()
        {
            try
            {
                File.Delete(fileName + "\\genotypes.inp");
                File.Delete(fileName + "\\haplotypes.out");
                File.Delete(fileName + "\\haplotypes.out_freqs");
                File.Delete(fileName + "\\haplotypes.out_hbg");
                File.Delete(fileName + "\\haplotypes.out_monitor");
                File.Delete(fileName + "\\haplotypes.out_pairs");
                File.Delete(fileName + "\\haplotypes.out_probs");
                File.Delete(fileName + "\\haplotypes.out_recom");
                return true;
            }
            catch (IOException e)
            {
                Console.WriteLine(Properties.Resources.PHASE_IMPROPERLY_TERMINATED + e);
                return false;
            }
            
        }

        #endregion

        #region good or bad

        /// <summary>
        /// Given a populated and calculated metric, determine based on handler's filter criteria whether
        /// that cluster is good or bad. Returns true for good, false for bad
        /// </summary>
        private bool GoodOrBad(ClusterMetric tempMetric)
        {
            if (tempMetric.Dirt > dirtCutoff
                    || tempMetric.AlignmentQuality < alignQualCutoff
                    || tempMetric.ReadQuality < readQualCutoff
                    || tempMetric.PopulationPercentage < populationPercentageCutoff 
                    || tempMetric.PloidyDisagreement > ploidyDisagreementCutoff)
            {
                tempMetric.Good = false;
            }
            else
            {
                tempMetric.Good = true;
            }
            return tempMetric.Good;
        }

        #endregion

        #region output files

        /// <summary>
        /// Add a sequence to the filtered output file header
        /// </summary>
        private void AddToHeader(SAMAlignedSequence seq)
        {
            newHeader.ReferenceSequences.Add(new ReferenceSequenceInfo(seq.RName, GetSequence(seq).Length)); 

            // for each good cluster
            SAMRecordField sq = new SAMRecordField("SQ");
            sq.Tags.Add(new SAMRecordFieldTag("SN", seq.RName));
            sq.Tags.Add(new SAMRecordFieldTag("LN", GetSequence(seq).Length.ToString(ci))); 
            newHeader.RecordFields.Add(sq);
        }

        private static string GetSequence(SAMAlignedSequence seq)
        {
            String seqStr = seq.QuerySequence.ToString();
            return Regex.Split(seqStr, "\r\n")[0]; 
        }

        /// <summary>
        /// Write metrics to all per-cluster metric files
        /// </summary>
        private void WriteToMetricOutputFiles(ClusterMetric metric, bool isGood)
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

        /// <summary>
        /// Write sequences from output queue into filtered bam output file
        /// </summary>
        /// <returns></returns>
        private System.Delegate WriteToBam()
        {
            while (bamOutputQueue.Count > 0)
            {
                Collection<SAMAlignedSequence> seqs = bamOutputQueue.Dequeue();
                foreach (SAMAlignedSequence seq in seqs)
                {
                    bamFormatter.WriteAlignedSequence(header, seq, bamStream);
                }
                seqs.Clear();
                seqs = null;
            }
            // signal to next thread runner that it can now process sequences from the queue
            canWriteToBam = true;
            return null;
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

        private void DeleteFilteredMetricFile()
        {
            formatterFilteredFile.Close();
            formatterFilteredFile = null;
            File.Delete(fileName + "\\filtered.metr");
            File.Delete(fileName + "\\filtered.bam");
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
                bamOutputQueue = new Queue<Collection<SAMAlignedSequence>>();

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
        /// Concatenate together the header and sequence bam output files
        /// </summary>
        private void concatFiles()
        {
            if (bamStream != null)
            {
                bamStream.Close();
                bamStream = null;
            }
            
            Console.WriteLine(Properties.Resources.CONSTRUCTING_OUTPUT_FILE);
            ProcessStartInfo startInfo = new ProcessStartInfo();

            startInfo.CreateNoWindow = true;
            startInfo.FileName = "cmd.exe";
            startInfo.Arguments = "/C copy /b " + fileName + "\\header.bam+" + fileName + "\\sequences.bam " + fileName + "\\filtered.bam /y";
            Console.WriteLine(Properties.Resources.COMMAND_LINE_INSTRUCTION + startInfo.Arguments);

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
                File.Delete(fileName + "\\header.bam");
                File.Delete(fileName + "\\sequences.bam");
            }
            
            bamFilesMerged = true;
            Console.WriteLine(Properties.Resources.CONSTRUCTING_OUTPUT_FILE_FINISHED);
        }

        #endregion

        #region summary and stats data

        private static void AddFxAverages(Collection<double> averagesList, Collection<double> frequencies)
        {
            while (averagesList.Count < frequencies.Count)
            {
                averagesList.Add(0);
            }
            int i = 0;
            foreach(double d in frequencies)
            {
                averagesList[i] += d;
                i++;
            }
        }

        /// <summary>
        /// Get various data arrays from metric, representing summary details for all clusters so far
        /// </summary>
        private void CreateSummaryArrays(ClusterMetric metric, bool isGood) {
                    
            clustSeqFrequencies.Add(metric.ClusterSequenceFrequencies);
            if(isGood) { clustSeqFrequenciesGood.Add(metric.ClusterSequenceFrequencies); }

            AddFxAverages(clusterSequenceFrequenciesOverview, metric.ClusterSequenceFrequencies);
            if (isGood) { AddFxAverages(clusterSequenceFrequenciesOverviewGood, metric.ClusterSequenceFrequencies); }
                    
            SetDictValueCounts(graphDataAllReads, metric.CountAll);
            if (isGood) { SetDictValueCounts(graphDataAllReadsGood, metric.CountAll); }
                    
            SetDictValueCounts(graphDataDistinctReads, metric.CountDistinct);
            if (isGood) { SetDictValueCounts(graphDataDistinctReadsGood, metric.CountDistinct); }
                    
            SetDictValueCounts(graphDataIndividualsCounts, metric.CountSamples);
            if(isGood) { SetDictValueCounts(graphDataIndividualsCountsGood, metric.CountSamples); }
                    
            int key = (int)Math.Round(metric.SampleReadCountsDistinct.Average(), 1);
            SetDictValueCounts(graphDataIndividualsDistinctReadcounts, key);
            if (isGood) { SetDictValueCounts(graphDataIndividualsDistinctReadcountsGood, key); }
                    
            key = (int)Math.Round(metric.SampleReadCountsAll.Average(), 1);
            SetDictValueCounts(graphDataIndividualsTotalReadcounts, key);
            if (isGood) { SetDictValueCounts(graphDataIndividualsTotalReadcountsGood, key); }
        }

        /// <summary>
        /// Set or update overview stats
        /// (Totals are used to enable easy obtaining of average without iterating through lists to count each time)
        /// </summary>
        private void SetOverviewStats(ClusterMetric metric, bool isGood)
        {
            // greatest number of samples found so far in any one cluster
            maxSampleCount = (metric.CountSamples > maxSampleCount) ? metric.CountSamples : maxSampleCount;

            readCountTotal += metric.CountAll;
            readCountDistinctTotal += metric.CountDistinct;
            if(isGood)
            {
                readCountGood += metric.CountAll;
                readCountDistinctGood += metric.CountDistinct;
            }

            // maximum quality value found so far in any one cluster
            maxMapQuality = (metric.AlignmentQuality > maxMapQuality) ? metric.AlignmentQuality : maxMapQuality;
            maxReadQuality = (metric.ReadQuality > maxReadQuality) ? metric.ReadQuality : maxReadQuality;

            // set totals for all clusters
            totalDirt += metric.Dirt;
            totalMapQ += metric.AlignmentQuality;
            totalReadQ += metric.ReadQuality;

            // set totals for good clusters
            totalDirtGood = (isGood) ? totalDirtGood + metric.Dirt : totalDirtGood;
            totalMapQGood += (isGood) ? metric.AlignmentQuality : 0;
            totalReadQGood += (isGood) ? metric.ReadQuality : 0;

            // running averages for all clusters
            averageDirt = Math.Round(totalDirt / (double)numberClustersParsed, 2);
            averageMapQ = Math.Round(totalMapQ / (double)numberClustersParsed, 2);
            averageReadQ = Math.Round(totalReadQ / (double)numberClustersParsed, 2);

            // running averages for good clusters
            averageDirtGood = (goodCount > 0) ? Math.Round(totalDirtGood / (double)goodCount, 2) : 0;
            averageMapQGood = (goodCount > 0) ? Math.Round(totalMapQGood / (double)goodCount, 2) : 0;
            averageReadQGood = (goodCount > 0) ? Math.Round(totalReadQGood / (double)goodCount, 2) : 0;
        }

        /// <summary>
        /// For the given key-count dictionary, if it has key, increment the value to reflect that it contains
        /// another instance. If it does not yet have that key, add a first instance with that key reference
        /// Dict values are rounded to the nearest whole number (int)
        /// </summary>
        private static void SetDictValueCounts(Dictionary<int, int> dict, int key)
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

        /// <summary>
        /// Outputs data to two files. Data represents read count (total, or distinct) and the number of clusters that have that count.
        /// These data files give us the actual, unnormalised size of each cluster (independent of which individuals the reads come from)
        /// The data could be used to determien if each cluster is above/below average size and/or the degree of distribution/variation
        /// in size between clusters
        /// 
        /// all_read_counts.tsv
        ///     total_num_reads | num_clusters
        /// 
        /// distinct_read_counts.tsv
        ///     num_distinct_reads | num_clusters
        /// </summary>
        private void PrintClustSummary_ReadCounts(Dictionary<int, int> allReads,
            Dictionary<int, int> distinctReads, Dictionary<int, int> individualsTotalReadcounts, 
            Dictionary<int, int> individualsDistinctReadcounts, string subdirectory)
        {
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "\\" + subdirectory + "\\" + "read_counts_all.tsv"))
            {
                file.WriteLine("#Read-per-cluster distribution, where read count includes all read instances including duplicates");
                file.WriteLine("total_num_reads\tnum_clusters");
                foreach (KeyValuePair<int, int> dat in allReads)
                {
                    file.WriteLine(dat.Key + "\t" + dat.Value);
                }
            }

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "\\" + subdirectory + "\\" + "read_counts_distinct.tsv"))
            {
                file.WriteLine("#Read-per-cluster distribution, where read count is distinct reads (duplicate reads will be counted as the one)");
                file.WriteLine("num_distinct_reads\tnum_clusters");
                foreach (KeyValuePair<int, int> dat in distinctReads)
                {
                    file.WriteLine(dat.Key + "\t" + dat.Value);
                }
            }

            // The below two files are the same as the above two, but read counts are the average per-individual count not
            // the total count in the cluster

            //(this report answers: on average, how many total read counts do the individuals in the cluster have?)
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "\\" + subdirectory + "\\" + "by_indiv_read_counts.tsv"))
            {
                file.WriteLine("#Read-count-per-individual distribution for each cluster, where read count includes all read instances including duplicates");
                file.WriteLine("total_num_reads\tnum_clusters");
                foreach (KeyValuePair<int, int> dat in individualsTotalReadcounts)
                {
                    file.WriteLine(dat.Key + "\t" + dat.Value);
                }
            }

            // average distinct readcount by individual 
            //(this report answers: on average, how many distinct read counts do the individuals in the cluster have)
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "\\" + subdirectory + "\\" + "by_indiv_read_counts_distinct.tsv"))
            {
                file.WriteLine("#Read-count-per-individual distribution for each cluster, where read count is distinct reads (duplicate reads will be counted as the one)");
                file.WriteLine("num_distinct_reads\tnum_clusters");
                foreach (KeyValuePair<int, int> dat in individualsDistinctReadcounts)
                {
                    file.WriteLine(dat.Key + "\t" + dat.Value);
                }
            }

        }

        /// <summary>
        /// The read-frequency distribution per-cluster (where frequency was calculated on a per-individual basis). This
        /// is the relative distribution of each sequence and the percentage which are the top-1, top-2,... top-n reads for
        /// each individual in the cluster
        /// Data is normalised and displayed as percentages from smaller to larger (e.g. piechart data)
        /// </summary>
        private void PrintClustSummary_FrequencyDistributions(List<Collection<double>> sequenceFrequencies, 
            string subdirectory)
        {
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "\\" + subdirectory + "\\" + "frequencies.tsv"))
            {
                file.WriteLine("#Per-cluster relative distribution of each sequence, from top-1, top-2, ...top-n");
                file.WriteLine("frequencies");
                foreach (Collection<double> dat in sequenceFrequencies)
                {
                    file.WriteLine(String.Join("\t", dat));
                }
            }
        }

        /// <summary>
        /// Per-cluster distribution of individuals (i.e. how many individuals are represented in each cluster, grouped by number of individuals)
        /// </summary>
        private void PrintClustSummary_IndivCounts(Dictionary<int, int> individualsCounts, string subdirectory)
        {
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "\\" + subdirectory + "\\" + "individuals_distribution.tsv"))
            {
                file.WriteLine("#Per-cluster distribution of individuals");
                file.WriteLine("num_individuals\tnum_clusters");
                foreach (KeyValuePair<int, int> dat in individualsCounts)
                {
                    file.WriteLine(dat.Key + "\t" + dat.Value);
                }
            }
        }

        // writes to a separate file for each significant all-clusters metric
        // also writes bam output file header (if required)
        private void PrintFinalCluserMetrics()
        {
            if (!finished) // this flag exists to make sure this method only executes once (unless reset)
            {

                // all sequences
                if(WriteOverviewMetricOriginal)
                {
                    if (!Directory.Exists(fileName + "\\" + OUTPUT_DIRECTORY))
                    {
                        Directory.CreateDirectory(fileName + "\\" + OUTPUT_DIRECTORY);
                    }
                    PrintClustSummary_ReadCounts(graphDataAllReads, graphDataDistinctReads, graphDataIndividualsTotalReadcounts, graphDataIndividualsDistinctReadcounts, OUTPUT_DIRECTORY);
                    PrintClustSummary_IndivCounts(graphDataIndividualsCounts, OUTPUT_DIRECTORY);
                    PrintClustSummary_FrequencyDistributions(clustSeqFrequencies, OUTPUT_DIRECTORY);
                }

                // filtered sequences
                if (WriteOverviewMetricFiltered && goodCount > 0)
                {
                    if (!Directory.Exists(fileName + "\\" + FILTERED_OUTPUT_DIRECTORY))
                    {
                        Directory.CreateDirectory(fileName + "\\" + FILTERED_OUTPUT_DIRECTORY);
                    }
                    PrintClustSummary_ReadCounts(graphDataAllReadsGood, graphDataDistinctReadsGood, 
                        graphDataIndividualsTotalReadcountsGood, graphDataIndividualsDistinctReadcountsGood, FILTERED_OUTPUT_DIRECTORY);
                    PrintClustSummary_IndivCounts(graphDataIndividualsCountsGood, FILTERED_OUTPUT_DIRECTORY);
                    PrintClustSummary_FrequencyDistributions(clustSeqFrequenciesGood, FILTERED_OUTPUT_DIRECTORY);
                }

                Console.WriteLine(Properties.Resources.FINISHED_WRITING_METRIC);
                finished = true;
            }

            while (!canWriteToBam)
            {
                Thread.Sleep(5000); // sleep 5 seconds, background thread is currently writing to bam file
            }
            
            // If bam output has been requested, and now that we know the file can be written to, continue
            if(writeToFilteredBam)
            {
                canWriteToBam = false;
                string headerFile = fileName + "\\header.bam";
                File.Delete(headerFile);    
                using (Stream writer = File.OpenWrite(headerFile))
                {
                    new BAMFormatter().WriteHeader(newHeader, writer);
                }
                canWriteToBam = true;
            }
        }

        #endregion

        #endregion


    }
}
