using Bio.Algorithms.Alignment;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace Bio.Algorithms.Metric
{
    /// <summary>
    /// An ClusterMetricHandler receives a set of SAMAlignedSequences (as a list or one by one)
    /// and passes them as "clusters" to a ClusterMetric, which performs a number of calculations
    /// on them to calculate the accuracy of the cluster.
    /// </summary>
    public class ClusterMetricHandler : IMetricHandler
    {
        #region Private Fields

        /// <summary>
        /// Name of file to which metric data will be written.
        /// </summary>
        private string fileName;

        /// <summary>
        /// Id of the cluster to which current sequences belong
        /// </summary>
        private string currCluster;

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
        
        private ClusterMetric metric = null;

        /// <summary>
        /// Used to write read data back out to SAM or BAM file
        /// </summary>
        private TextWriter samWriter = null;

        /// <summary>
        /// Indicates whether input is to be written back out to a SAM or BAM file
        /// </summary>
        private bool writeToFile = false;

        #endregion

        #region Constructors

        /// <summary>
        /// The default constructor.
        /// </summary>
        public ClusterMetricHandler()
        {
            //throw new NotImplementedException();
        }

        /// <summary>
        /// Non-default constructor, used to set the file name
        /// </summary>
        /// <param name="clusterFileName">Name of the file metric data is to be written to.</param>
        public ClusterMetricHandler(string clusterFileName)
            : this()
        {
            FileName = clusterFileName;
        }

        /// <summary>
        /// Non-default constructor, used to set the file name and also add a list of sequences
        /// </summary>
        /// <param name="clusterFileName">Name of the file metric data is to be written to.</param>
        /// <param name="sequences">List of sequences.</param>
        public ClusterMetricHandler(string clusterFileName, List<SAMAlignedSequence> sequences)
            : this(clusterFileName)
        {
            AddAll(sequences);
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
        public void AddAll(List<SAMAlignedSequence> sequences)
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
                if (metric == null)
                {
                    metric = new ClusterMetric();
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

                metric.Calculate(sequences);
                Console.WriteLine(metric.ToFileString()); // todo fixme this line of course should be removed

                formatter.Write(metric);

                if (writeToFile) // todo if the metric determines that we want to keep this read for downstream analysis
                {
                    foreach (IAlignedSequence seq in sequences)
                    {
                        SAMFormatter.WriteSAMAlignedSequence(seq, samWriter);
                    }
                }
                metric.Reset();
            }
        }

        /// <summary>
        /// Process any/all remaining sequences.
        /// </summary>
        public void FlushSequences()
        {
            ProcessSequences();
            // Any other cleanup required? todo
        }

        /// <summary>
        /// Disposes all formatters
        /// </summary>
        public void Dispose()
        {
            FlushSequences();
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
