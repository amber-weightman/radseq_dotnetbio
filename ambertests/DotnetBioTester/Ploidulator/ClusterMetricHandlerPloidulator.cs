using Bio.Algorithms.Alignment;
using Bio.Algorithms.Metric;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Text;
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
        /// GUI dispatcher, used to queue jobs on GUI thread
        /// </summary>
        private Dispatcher dispatcher;

        /// <summary>
        /// Panel in which to display visualisation per cluster
        /// </summary>
        private WrapPanel panel;

        private List<ClusterMetricPloidulator> clusterBucket;

        /// <summary>
        /// Id of the cluster to which current sequences belong
        /// </summary>
        private string currCluster;

        /// <summary>
        /// Number of clusters which have been received
        /// </summary>
        private int clusterCount = 0;

        /// <summary>
        /// Sequences stored, not yet processed.
        /// </summary>
        private List<SAMAlignedSequence> sequences = null;

        /// <summary>
        /// Every sequence, for if we want to store in mem
        /// </summary>
        private List<SAMAlignedSequence> allSequences = null;
        // curr sequences
        //private SequenceAlignmentMap seqMap = null;

        /// <summary>
        /// Formatter used to write metric data to file
        /// </summary>
        private MetricFormatter formatterOriginalFile = null;
        private MetricFormatter formatterFilteredFile = null;

        /// <summary>
        /// ClusterMetric calculates metric values for each SAMAlignedSequence list
        /// </summary>

        private ClusterMetricPloidulator metric = null;

        /// <summary>
        /// Used to write read data back out to SAM or BAM file
        /// </summary>
        private TextWriter samWriter = null;
        private Stream bamWriter = null;

        private BAMFormatter bamFormatter = null;

        private SAMAlignmentHeader header = null;

        /// <summary>
        /// Indicates whether input is to be written back out to a SAM or BAM file
        /// </summary>
        private bool writeToFilteredBam = true;
        private bool writeClusterMetricOriginal = true;
        private bool writeClusterMetricFiltered = true;
        private bool writeOverviewMetricOriginal = true;
        private bool writeOverviewMetricFiltered = true;
        public bool WriteToFilteredBam { get { return writeToFilteredBam; } set { writeToFilteredBam = value; } }
        public bool WriteClusterMetricOriginal { get { return writeClusterMetricOriginal; } set { writeClusterMetricOriginal = value; } }
        public bool WriteClusterMetricFiltered { get { return writeClusterMetricFiltered; } set { writeClusterMetricFiltered = value; } }
        public bool WriteOverviewMetricOriginal { get { return writeOverviewMetricOriginal; } set { writeOverviewMetricOriginal = value; } }
        public bool WriteOverviewMetricFiltered { get { return writeOverviewMetricFiltered; } set { writeOverviewMetricFiltered = value; } }


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

        double dirtCutoff = 1;
        double alignQualCutoff = 0;
        double readQualCutoff = 0;
        double popPercent = 0;
        int numSamples = 0;
        int numClustersToParse = int.MaxValue;

        private int numClustersParsed = 0;
        private int goodCount = 0;
        private int maxSampleCount = 0;
        private double maxMapQ = 0;
        private double maxReadQ = 0;
        private double totalDirt = 0;
        private double totalMapQ = 0;
        private double totalReadQ = 0;
        private double totalDirtGood = 0;
        private double totalMapQGood = 0;
        private double totalReadQGood = 0;

        private double averageDirt = 0;
        private double averageMapQ = 0;
        private double averageReadQ = 0;
        private double averageDirtGood = 0;
        private double averageMapQGood = 0;
        private double averageReadQGood = 0;

        /// <summary>
        /// precondition: metric.Good must have been set
        /// </summary>
        /// <param name="metric"></param>
        private void SetStats(ClusterMetricPloidulator metric)
        {
            // greatest number of samples found so far in any one cluster
            maxSampleCount = (metric.CountSamples > maxSampleCount) ? metric.CountSamples : maxSampleCount;

            // maximum quality value found so far in any one cluster
            maxMapQ = (metric.AlignmentQualities[0] > maxMapQ) ? metric.AlignmentQualities[0] : maxMapQ;
            maxReadQ = (metric.ReadQualities[0] > maxReadQ) ? metric.ReadQualities[0] : maxReadQ;

            // totals are used to enable easy obtaining of average without iterating through lists to count each time
            totalDirt += metric.Dirt;
            totalMapQ += metric.AlignmentQualities[0];
            totalReadQ += metric.ReadQualities[0];

            totalDirtGood = (metric.Good) ? totalDirtGood + metric.Dirt : totalDirtGood;
            
            totalMapQGood += (metric.Good) ? metric.AlignmentQualities[0] : 0;
            totalReadQGood += (metric.Good) ? metric.ReadQualities[0] : 0;

            // running-averages
            averageDirt = Math.Round(totalDirt / (double)numClustersParsed, 2);
            averageMapQ = Math.Round(totalMapQ / (double)numClustersParsed, 2);
            averageReadQ = Math.Round(totalReadQ / (double)numClustersParsed, 2);

            averageDirtGood = (goodCount > 0) ? Math.Round(totalDirtGood / (double)goodCount, 2) : 0;
            averageMapQGood = (goodCount > 0) ? Math.Round(totalMapQGood / (double)goodCount, 2) : 0;
            averageReadQGood = (goodCount > 0) ? Math.Round(totalReadQGood / (double)goodCount, 2) : 0;
        }


        #endregion

        #region Constructors

        /// <summary>
        /// The default constructor.
        /// </summary>
        public ClusterMetricHandlerPloidulator()
        {
            //throw new NotImplementedException();
        }

        /// <summary>
        /// Non-default constructor, used to set the file name
        /// </summary>
        /// <param name="clusterFileName">Name of the file metric data is to be written to.</param>
        public ClusterMetricHandlerPloidulator(string clusterFileName)
            : this()
        {
            FileName = clusterFileName;
        }

        /// <summary>
        /// Non-default constructor, used to set the file name
        /// </summary>
        /// <param name="clusterFileName">Name of the file metric data is to be written to.</param>
        public ClusterMetricHandlerPloidulator(string clusterFileName, int ploidy, 
            double dirtCutoff, double alignQualCutoff, double readQualCutoff, double popPercent, int numSamples, bool? outputToFile,
            Dispatcher d, WrapPanel p)
            : this()
        {
            FileName = clusterFileName;
            expectedPloidy = ploidy;
            dispatcher = d;
            panel = p;

            this.dirtCutoff = dirtCutoff;
            this.alignQualCutoff = alignQualCutoff;
            this.readQualCutoff = readQualCutoff;
            this.popPercent = popPercent;
            this.numSamples = numSamples;
            this.writeToFilteredBam = (outputToFile.HasValue) ? outputToFile.Value : false;

        }



        /// <summary>
        /// Non-default constructor, used to set the file name and also add a list of sequences
        /// </summary>
        /// <param name="clusterFileName">Name of the file metric data is to be written to.</param>
        /// <param name="sequences">List of sequences.</param>
        public ClusterMetricHandlerPloidulator(string clusterFileName, Collection<SAMAlignedSequence> sequences)
            : this(clusterFileName)
        {
            AddRange(sequences);
        }

        #endregion

        #region Properties

        /// <summary>
        /// Name of file to which metric data will be written.
        /// </summary>
        public string FileName
        {
            get { return fileName; }
            set { this.fileName = value; }
        }

        /// <summary>
        /// Id of the cluster to which current sequences belong.
        /// </summary>
        public string CurrCluster
        {
            get { return currCluster; }
        }

        /// <summary>
        /// number of clusters which have been parsed
        /// </summary>
        public int ClusterCount
        {
            get { return clusterCount; }
        }


        /// <summary>
        /// Sequences stored, not yet processed.
        /// </summary>
        public List<SAMAlignedSequence> Sequences
        {
            get { return sequences; }
        }

        /// <summary>
        /// Average frequencies for each sequence, calculated per sample per cluster and averaged
        /// out to cluster
        /// </summary>
        public List<double[]> ClustSeqFrequencies
        {
            get { return clustSeqFrequencies; }
        }


        public SAMAlignmentHeader Header
        {
            get { return header; }
            set { header = value; }
        }


        

        public int GoodCount { get { return goodCount; } }
        public int MaxSampleCount { get { return maxSampleCount; } }
        public double MaxMapQ { get { return maxMapQ; } }
        public double MaxReadQ { get { return maxReadQ; } }
        public double AverageDirt { get { return averageDirt; } }
        public double AverageMapQ { get { return averageMapQ; } }
        public double AverageReadQ { get { return averageReadQ; } }

        public double AverageDirtGood { get { return averageDirtGood; } }
        public double AverageMapQGood { get { return averageMapQGood; } }
        public double AverageReadQGood { get { return averageReadQGood; } }

        private bool aborted = false;

        public void Abort()
        {
            aborted = true;
        }

        #endregion

        #region Public Methods

        // this needs a couple more params 
        // -- one to specify whether to write to bam or sam
        // -- another param which is the filter criteria (filter reads out based on cluster dirt at what cutoff)
        /// <param name="writeToFile">Whether the SAMAlignedSequences, filtered by the ClusterMetric, should
        /// be written back out to a new SAM/BAM file.</param>
        public void WriteToFile(bool writeToFile)
        {
            this.writeToFilteredBam = writeToFile;
        }

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
            if (!finished)
            {
                string thisSeqCluster = sequence.RName; // Cluster the sequence we just added belongs to

                if (currCluster == null)
                {
                    currCluster = thisSeqCluster;
                }
                else if (!currCluster.Equals(thisSeqCluster)) // This sequence belongs to a different cluster from the ones currently stored in sequences
                {
                    ++numClustersParsed; // mark off another cluster
                    ProcessSequences();
                    
                    currCluster = thisSeqCluster;
                    sequences = new List<SAMAlignedSequence>();
                }

                if (sequences == null)
                {
                    sequences = new List<SAMAlignedSequence>();
                }
                sequences.Add(sequence);
                return true;
            }
            // we should be finished but we are still outputting to the bam file
            else if(!canWriteToBam)
            {
                Thread.Sleep(5000); // sleep 5 seconds
                return true;
            }
            return false;
        }

        

        /// <summary>
        /// Calculate metric for a single cluser of sequences (all stored sequences), 
        /// and write metric data to file.
        /// </summary>
        public void ProcessSequences()
        {
            if (sequences != null && sequences.Count > 0)
            {
                //if ((metric == null || metric.Id != "142919") && (clusterCount < numClustersToParse)) 
                //{
                    clusterCount++;
                    Console.Write("a");
                
                    /*if (metric == null)
                    {
                        metric = new ClusterMetricPloidulator(expectedPloidy, numSamples, dispatcher, panel);
                    }*/

                    ClusterMetricPloidulator tempMetric = new ClusterMetricPloidulator(expectedPloidy, numSamples, dispatcher, panel);

                    if (writeClusterMetricOriginal && formatterOriginalFile == null)
                    {
                        formatterOriginalFile = new MetricFormatter(FileName + "_orig.metr");
                    }
                    if (writeClusterMetricFiltered && formatterFilteredFile == null)
                    {
                        formatterFilteredFile = new MetricFormatter(FileName + "_filtered.metr");
                    }
                    Console.Write("b");
                    if (writeToFilteredBam && bamWriter == null)
                    {
                        bamWriter = File.Create(FileName + "_filtered.bam");
                        bamFormatter = new BAMFormatter();   
                        bamFormatter.WriteHeader(header, bamWriter);
                        bamOutputQueue = new Queue<List<SAMAlignedSequence>>();

                    }
                    Console.Write("c");
                    // calculate metric
                    //metric.Calculate(sequences);

                    


                    tempMetric.Calculate(sequences);

                    Console.Write("d");
                    // read and/or store some values from the metric
                    clustSeqFrequencies.Add(tempMetric.ClustSeqFrequencies);
                    Console.Write("e");
                    SetThisDictValueCounts(graphDataAllReads, tempMetric.CountAll);
                    Console.Write("f");
                    SetThisDictValueCounts(graphDataDistinctReads, tempMetric.CountDistinct);
                    Console.Write("g");
                    SetThisDictValueCounts(graphDataIndividualsCounts, tempMetric.CountSamples);
                    Console.Write("h");
                    int key = (int)Math.Round(tempMetric.SampleReadCountsDistinct.Average(), 1);
                    SetThisDictValueCounts(graphDataIndividualsDistinctReadcounts, key);
                    Console.Write("i");
                    key = (int)Math.Round(tempMetric.SampleReadCountsAll.Average(), 1);
                    SetThisDictValueCounts(graphDataIndividualsTotalReadcounts, key);
                    Console.Write("j");
                    bool isOk = true;
                    Console.Write("k");
                    Console.WriteLine(tempMetric.PopulationPercentage + "********************");
                    if (tempMetric.Dirt > dirtCutoff 
                        || tempMetric.AlignmentQualities[0] < alignQualCutoff 
                        || tempMetric.ReadQualities[0] < readQualCutoff
                        || tempMetric.PopulationPercentage < popPercent)
                    {
                        isOk = false;
                    }
                    if (isOk)
                    {
                        ++goodCount;
                        Console.WriteLine("GOOD "+goodCount);
                    }
                    else
                    {
                        Console.WriteLine("BAD");
                    }
                    tempMetric.Good = isOk;
                    SetStats(tempMetric);

                    Console.Write("<l>");
                    // write metric out to metric file
                    if (writeClusterMetricOriginal)
                    {
                        formatterOriginalFile.Write(tempMetric);
                    }
                    if (writeClusterMetricFiltered && isOk)
                    {
                        formatterFilteredFile.Write(tempMetric);
                    }
                    
                    Console.Write("m");

                    
                    

                    // if cluster is good, write aligned reads out to another sam/bam file for downstream analysis
                    if (writeToFilteredBam && isOk)
                    {
                        Console.WriteLine("writing good cluster to file for sequences: " + sequences.Count);
                        while (bamOutputQueue.Count > 100)
                        {
                            // Gives the slow bamWriter a chance to catch up 
                            // fixme - some better way of checking this.
                            Thread.Sleep(20000); // sleep 20 seconds
                        }
                        bamOutputQueue.Enqueue(sequences);
                    } 
                    // we might be writing these to file later, so store them
                    // we might be changing the parameters for isOk later, so for now all is ok
                    else if(!writeToFilteredBam)
                    {
                        if(clusterBucket == null)
                        {
                            clusterBucket = new List<ClusterMetricPloidulator>();
                        }
                        StoreMetric(tempMetric);
                    }
                    if (writeToFilteredBam && canWriteToBam && bamOutputQueue.Count > 0)
                    {
                        canWriteToBam = false;
                        runner = new SequenceDelegate(WriteToBam);
                        runner.BeginInvoke(null, null);
                    }
                    
                    
                //}
                //else if (numClustersToParse != -1 || metric.Id == "142919") // for debugging pretend this is the end of the file if numClustersToParse is not -1
                //{
                  //  PrintAllCluserMetrics();
                //}
            }
            if(aborted)
            {
                // todo this might not quite tidy it up the way we want it to
                PrintAllCluserMetrics();
            }
        }
        private bool canWriteToBam = true;
        private bool finished = false;
        private void StoreMetric(ClusterMetricPloidulator tempMetric)
        {
            Console.Write("o");
            tempMetric.Shrink();
            Console.Write("p");
            clusterBucket.Add(tempMetric);
            Console.Write("q");
        }

        /// <summary>
        /// this is the poisson-graphing type integer dict
        /// all values are rounded to the nearest whole numbers
        /// </summary>
        /// <param name="dict"></param>
        /// <param name="key"></param>
        private void SetThisDictValueCounts(Dictionary<int, int> dict, int key)
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

        private void PrintAllCluserMetrics()
        {
            if (!finished)
            {


                // this assumes at the moment that no clusters have been filtered out


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
        }

        public delegate System.Delegate OneArgDelegate(double[] data);
        public delegate System.Delegate MetricDelegate(ClusterMetricPloidulator metric);
        public delegate System.Delegate SequenceDelegate();

        public SequenceDelegate runner;
        private Queue<List<SAMAlignedSequence>> bamOutputQueue;

        private int numChartsDisplayed = 0;
        private int maxNumCharts = 10;
        public int MaxNumCharts
        {
            get { return maxNumCharts; }
            set { this.maxNumCharts = value; }
        }
        private static int CHART_INCREMENT = 10; // syntax?


        public System.Delegate WriteToBam()
        {
            while (bamOutputQueue.Count > 0)
            {
                int i = 0;
                foreach (SAMAlignedSequence seq in bamOutputQueue.Dequeue())
                {
                    //Console.Write("-"+(i++)+"-");
                    bamFormatter.WriteAlignedSequence(header, seq, bamWriter);
                }
            }

            // signal to next thread runner that it can now process
            canWriteToBam = true;
            return null;
        }

        



        /// <summary>
        /// Process any/all remaining sequences.
        /// </summary>
        public void SetComplete()
        {
            ProcessSequences();

            PrintAllCluserMetrics();
            Console.WriteLine("im finished that, so what else do you want to do?");
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
            if (samWriter != null)
            {
                samWriter.Close();
            }
            if (bamWriter != null)
            {
                bamWriter.Close();
            }
            if (bamFormatter != null)
            {
                //bamFormatter.Close();
            }
        }

        #endregion

    }
}
