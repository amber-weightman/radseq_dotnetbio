using Bio.Algorithms.Alignment;
using Bio.Algorithms.Metric;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Text;
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

        /// <summary>
        /// Id of the cluster to which current sequences belong
        /// </summary>
        private string currCluster;

        /// <summary>
        /// Number of clusters which have been received
        /// </summary>
        private int clusterCount;

        /// <summary>
        /// Sequences stored, not yet processed.
        /// </summary>
        private List<SAMAlignedSequence> sequences = null;

        /// <summary>
        /// Formatter used to write metric data to file
        /// </summary>
        private MetricFormatter formatter = null;

        /// <summary>
        /// ClusterMetric calculates metric values for each SAMAlignedSequence list
        /// </summary>

        private ClusterMetricPloidulator metric = null;

        /// <summary>
        /// Used to write read data back out to SAM or BAM file
        /// </summary>
        private TextWriter samWriter = null;

        /// <summary>
        /// Indicates whether input is to be written back out to a SAM or BAM file
        /// </summary>
        private bool writeToFile = false;

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
        int alignQualCutoff = 0;
        string readQualCutoff = "na unused placeholder";
        int numSamples = 0;
        int numClustersToParse = int.MaxValue;


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
            double dirtCutoff, int alignQualCutoff, string readQualCutoff, int numSamples, int numClustersToParse,
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
            this.numSamples = numSamples;
            this.numClustersToParse = numClustersToParse;
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

        #endregion

        #region Public Methods

        // this needs a couple more params 
        // -- one to specify whether to write to bam or sam
        // -- another param which is the filter criteria (filter reads out based on cluster dirt at what cutoff)
        /// <param name="writeToFile">Whether the SAMAlignedSequences, filtered by the ClusterMetric, should
        /// be written back out to a new SAM/BAM file.</param>
        public void WriteToFile(bool writeToFile)
        {
            this.writeToFile = writeToFile;
        }

        /// <summary>
        /// Add a list of sequences. For each sequence, if the sequence belongs to the current cluster, 
        /// store it. If the sequence is part of a new cluster, process the current sequence cluster 
        /// then add the sequenceto a new cluster
        /// </summary>
        /// <param name="sequences">A list of sequences.</param>
        public void AddRange(IEnumerable<SAMAlignedSequence> sequences)
        {
            foreach (SAMAlignedSequence seq in sequences)
            {
                Add(seq);
            }
        }

        /// <summary>
        /// Add a sequence. If the sequence belongs to the current cluster, store it. If the sequence
        /// is part of a new cluster, process the current sequence cluster then add the sequence
        /// to a new cluster
        /// </summary>
        /// <param name="sequence">A sequence.</param>
        public void Add(SAMAlignedSequence sequence)
        {
            string thisSeqCluster = sequence.RName; // Cluster the sequence we just added belongs to
            if (currCluster == null)
            {
                currCluster = thisSeqCluster;
            }
            else if (!currCluster.Equals(thisSeqCluster)) // This sequence belongs to a different cluster from the ones currently stored in sequences
            {
                ProcessSequences();
                currCluster = thisSeqCluster;
                sequences = new List<SAMAlignedSequence>();
            }

            if (sequences == null)
            {
                sequences = new List<SAMAlignedSequence>();
            }
            sequences.Add(sequence);
        }

        /// <summary>
        /// Calculate metric for a single cluser of sequences (all stored sequences), 
        /// and write metric data to file.
        /// </summary>
        public void ProcessSequences()
        {
            if (sequences != null && sequences.Count > 0)
            {
                if (clusterCount < numClustersToParse || numClustersToParse == -1) 
                {
                    clusterCount++;

                
                    if (metric == null)
                    {
                        metric = new ClusterMetricPloidulator(expectedPloidy, numSamples, dispatcher, panel);
                    }
                    if (formatter == null)
                    {
                        formatter = new MetricFormatter(FileName + ".metr");
                    }
                    if (writeToFile && samWriter == null)
                    {
                        // todo fixme we might prefer this to be a bam file, but i'm using sam now so I can easily read it
                        samWriter = new StreamWriter(FileName + "_filtered.sam");
                        SAMFormatter.WriteHeader(new SAMAlignmentHeader(), samWriter); // todo what goes in the header?
                    }

                    // calculate metric
                    metric.Calculate(sequences);

                    // read and/or store some values from the metric
                    clustSeqFrequencies.Add(metric.ClustSeqFrequencies);

                    SetThisDictValueCounts(graphDataAllReads, metric.CountAll);
                    /*if (graphDataAllReads.ContainsKey(metric.CountAll))
                    {
                        ++graphDataAllReads[metric.CountAll];
                    }
                    else
                    {
                        graphDataAllReads[metric.CountAll] = 1;
                    }*/

                    SetThisDictValueCounts(graphDataDistinctReads, metric.CountDistinct);
                    /*if (graphDataDistinctReads.ContainsKey(metric.CountDistinct))
                    {
                        ++graphDataDistinctReads[metric.CountDistinct];
                    }
                    else
                    {
                        graphDataDistinctReads[metric.CountDistinct] = 1;
                    }*/

                    SetThisDictValueCounts(graphDataIndividualsCounts, metric.CountSamples);
                    /*if (graphDataIndividualsCounts.ContainsKey(metric.CountSamples))
                    {
                        ++graphDataIndividualsCounts[metric.CountSamples];
                    }
                    else
                    {
                        graphDataIndividualsCounts[metric.CountSamples] = 1;
                    }*/

                    int key = (int)Math.Round(metric.SampleReadCountsDistinct.Average(), 1);
                    SetThisDictValueCounts(graphDataIndividualsDistinctReadcounts, key);
                    /*if (graphDataIndividualsDistinctReadcounts.ContainsKey(key))
                    {
                        ++graphDataIndividualsDistinctReadcounts[key];
                    }
                    else
                    {
                        graphDataIndividualsDistinctReadcounts[key] = 1;
                    }*/
                    

                    key = (int)Math.Round(metric.SampleReadCountsAll.Average(), 1);
                    SetThisDictValueCounts(graphDataIndividualsTotalReadcounts, key);
                    /*if (graphDataIndividualsTotalReadcounts.ContainsKey(key))
                    {
                        ++graphDataIndividualsTotalReadcounts[key];
                    }
                    else
                    {
                        graphDataIndividualsTotalReadcounts[key] = 1;
                    }*/

                    bool isOk = true;
                    if (metric.Dirt > dirtCutoff)
                    {
                        isOk = false;
                    }

                    if (metric.AlignmentQual < alignQualCutoff)
                    {
                        isOk = false;
                    }
                    if(!isOk)
                    {
                        //Console.WriteLine("BAD");
                    }

                    // debugging purposes only
                    if (Convert.ToInt32(metric.Id) >= numClustersToParse) // only works when clust id is a num
                    {
                        // end iterating and calculate metrics across clusters
                    }

                    // write metric out to metric file
                    formatter.Write(metric);

                    // if cluster is good, write aligned reads out to another sam/bam file for downstream analysis
                    if (writeToFile && metric.Good) 
                    {
                        foreach (IAlignedSequence seq in sequences)
                        {
                            SAMFormatter.WriteSAMAlignedSequence(seq, samWriter);
                        }
                    }

                    // draw a chart for this cluster
                    /*for (; numChartsDisplayed < maxNumCharts && clustSeqFrequencies.Count > numChartsDisplayed; numChartsDisplayed++)
                    {
                        dispatcher.BeginInvoke(
                        System.Windows.Threading.DispatcherPriority.Normal,
                        new OneArgDelegate(DrawChart),
                        clustSeqFrequencies[numChartsDisplayed]);
                    }*/
  
                    metric.Reset();
                }
                else if (numClustersToParse != -1) // for debugging pretend this is the end of the file if numClustersToParse is not -1
                {
                    PrintAllCluserMetrics();

                }
            }
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
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_allreadcounts"))
            {
                file.WriteLine("num_reads\tnum_clusters");
                foreach (KeyValuePair<int, int> dat in graphDataAllReads)
                {
                    file.WriteLine(dat.Key + "\t" + dat.Value);
                }
            }

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_distinctreadcounts"))
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

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_numindividualscounts"))
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
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_byindividualsrtotaleadcounts"))
            {
                file.WriteLine("num_reads\tnum_clusters");
                foreach (KeyValuePair<int, int> dat in graphDataIndividualsTotalReadcounts)
                {
                    file.WriteLine(dat.Key + "\t" + dat.Value);
                }
            }

            // average distinct readcount by individual
            //graphDataIndividualsDistinctReadcounts
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_byindividualsrdistincteadcounts"))
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
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName + "_frequencies"))
            {
                file.WriteLine("frequencies");
                foreach (double[] dat in clustSeqFrequencies)
                {
                    file.WriteLine(String.Join("\t", dat));
                }
            }
            #endregion




            Console.WriteLine("finished");
        }

        public delegate System.Delegate OneArgDelegate(double[] data);
        private int numChartsDisplayed = 0;
        private int maxNumCharts = 10;
        public int MaxNumCharts
        {
            get { return maxNumCharts; }
            set { this.maxNumCharts = value; }
        }
        private static int CHART_INCREMENT = 10; // syntax?


        public System.Delegate DrawChart(double[] data)
        {
            SolidColorBrush myTransparentBrush = new SolidColorBrush();
            myTransparentBrush.Color = (Color)ColorConverter.ConvertFromString("Transparent");

            // Chart frame
            Chart c = new Chart();
            c.Name = "PieChart";
            c.VerticalAlignment = VerticalAlignment.Top;
            c.Margin = new System.Windows.Thickness(0, 0, 0, 0);
            c.Height = 200;
            c.Width = 200;
            c.Foreground = myTransparentBrush;
            c.BorderBrush = myTransparentBrush;

            // Chart data
            PieSeries s = new PieSeries();
            s.Title = "Frequencies";
            s.Margin = new System.Windows.Thickness(0, -40, 0, -20);
            c.Foreground = myTransparentBrush;
            c.Background = myTransparentBrush;
            c.BorderBrush = myTransparentBrush;
            s.IndependentValueBinding = new Binding("Key");
            s.DependentValueBinding = new Binding("Value");
            c.Series.Add(s);

            // Legend style
            Style style = new Style(typeof(Control));
            style.Setters.Add(new Setter(Chart.HeightProperty, new Binding("0")));
            style.Setters.Add(new Setter(Chart.WidthProperty, new Binding("0")));
            style.Setters.Add(new Setter(Chart.BorderBrushProperty, new Binding("Transparent")));
            style.Setters.Add(new Setter(Chart.BackgroundProperty, new Binding("Transparent")));
            c.LegendStyle = style;

            // Bind data
            KeyValuePair<string, double>[] dataSeries = new KeyValuePair<string, double>[data.Length];
            int i = 0;
            foreach (double d in data)
            {
                dataSeries[i] = new KeyValuePair<string, double>("Ploidy " + (++i), d);
            }
            ((PieSeries)c.Series[0]).ItemsSource = dataSeries;

            // Add chart
            panel.Children.Add(c);
            return null;
        }



        /// <summary>
        /// Process any/all remaining sequences.
        /// </summary>
        public void SetComplete()
        {
            ProcessSequences();
            PrintAllCluserMetrics();
            // Any other cleanup required that doesn't fit into Dispose()? todo
        }

        /// <summary>
        /// Disposes all formatters
        /// </summary>
        public void Dispose()
        {
            SetComplete();
            // todo fixme i have probably done this wrong. there will be some other checking and stuff
            // that i also need to do or other stuff I need to clean up
            if (formatter != null)
            {
                formatter.Close();
            }
            if (samWriter != null)
            {
                samWriter.Close();
            }
        }

        #endregion

    }
}
