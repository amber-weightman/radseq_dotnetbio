using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Bio.Algorithms.Metric
{
    /// <summary>
    /// An IMetric calculates various metric values from a list of SAMAlignedSequences, where each
    /// list of sequences belongs to a cluster. Clusters are determined based on sequence similarity, 
    /// e.g. as calculated by an MCL graph clustering algorithm, or any other clustering approach
    /// </summary>
    class ClusterMetric : IMetric
    {

        #region Private Fields

        /// <summary>
        /// A dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value
        /// </summary>
        private Dictionary<String, List<SAMAlignedSequence>> sequenceDict = null;

        /// <summary>
        /// A dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// </summary>
        private Dictionary<String, List<SAMAlignedSequence>> sampleDict = null;

        /// <summary>
        /// List of all sequences.
        /// </summary>
        private List<SAMAlignedSequence> sequences = null;

        #endregion

        #region Constructors

        /// <summary>
        /// The default constructor.
        /// </summary>
        public ClusterMetric()
        {
            //throw new NotImplementedException();
        }

        #endregion

        #region Properties

        /// <summary>
        /// Cluster ID (unique ID of reference sequence against which all sequences in the cluster
        /// are aligned)
        /// </summary>
        public string Id
        {
            // All sequences are assumed to be from the same cluster/chromosome. todo fixme this is potentially
            // not good if we want to be able to re-cluster sequences?
            get 
            {
                string rname = sequences[0].RName;
                /*
                byte[] bytes = new byte[rname.Length * sizeof(char)];
                System.Buffer.BlockCopy(rname.ToCharArray(), 0, bytes, 0, bytes.Length);

                if (bytes[0x4] == 0)
                {
                    // so what do I actually do?
                }*/

                return (sequences != null) ? rname : null; 
            }

            /*
             * One final thought: just taking RNAMEs doesn't always exactly indicate reference 
             * assignment/cluster membership.  In the past, I've also filtered on flag 4, from this 
             * paragraph in the sam specification:

                "Bit 0x4 is the only reliable place to tell whether the read is unmapped. If 0x4 is set, no
                assumptions can be made about RNAME, POS, CIGAR, MAPQ, bits 0x2, 0x10, 0x100 and
                0x800, and the bit 0x20 of the previous read in the template."
             * */

        }

        /// <summary>
        /// Total number of all sequences in the cluster
        /// </summary>
        public int CountAll
        {
            get { return (sequences != null) ? sequences.Count : 0; }
        }

        /// <summary>
        /// Total number of distinct sequences in the cluster
        /// </summary>
        public int CountDistinct
        {

            get { return (sequenceDict != null) ? sequenceDict.Count : 0; }
        }

        /// <summary>
        /// Total number of samples represented in the cluster
        /// </summary>
        public int CountSamples
        {
            get { return (sampleDict != null) ? sampleDict.Count : 0; }
        }

        /// <summary>
        /// A dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value
        /// </summary>
        public Dictionary<String, List<SAMAlignedSequence>> SequenceDict
        {
            get { return sequenceDict; }
        }

        /// <summary>
        /// A dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// </summary>
        public Dictionary<String, List<SAMAlignedSequence>> SampleDict
        {
            get { return sampleDict; }
        }

        /// <summary>
        /// List of all sequences.
        /// </summary>
        public List<SAMAlignedSequence> Sequences
        {
            get { return sequences; }
        }

        /// <summary>
        /// The total number of sequences which are within ploidy (i.e. the top two represented query
        /// sequences)
        /// </summary>
        public int CountInPloidy
        {
            get
            {
                int count = 0, seq0 = 0, seq1 = 0;
                // todo fixme - i am going to say this could be what is not working. 
                // the dict is not reliably sorted?
                var sortedSequenceDict = (from sequence in sequenceDict
                                          orderby sequence.Value.Count descending
                                          select sequence)
                    .ToDictionary(pair => pair.Key, pair => pair.Value);

                foreach (KeyValuePair<String, List<SAMAlignedSequence>> seqs in sortedSequenceDict)
                {
                    if (count == 0)
                    {
                        seq0 = seqs.Value.Count;
                        ++count;
                    }
                    else if (count == 1)
                    {
                        seq1 = seqs.Value.Count;
                        ++count;
                    }
                    else
                    {
                        break;
                    }
                }
                return seq0 + seq1;
            }
        }

        /// <summary>
        /// The total number of sequences which are beyond ploidy (i.e. total count excluding the 
        /// top two represented query sequences)
        /// </summary>
        public Double BeyondPloidy
        {
            get { return (Double)(CountAll - CountInPloidy) / (Double)CountAll; }
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// String of tab-separated values for writing to file by MetricFormatter.
        /// todo fixme we want the MetricFormater to decide how to write each line, but
        /// for now until I know what values are being written I will let the Metric handle this.
        /// </summary>
        public string ToFileString()
        {
            return Id + "\t" +
                CountDistinct + "\t" +
                CountSamples + "\t" +
                BeyondPloidy;
        }

        /// <summary>
        /// Calculate metric values from the given list of sequences.
        /// </summary>
        public void Calculate(List<SAMAlignedSequence> sequences)
        {
            this.sequences = sequences;
            SetSequenceDict();
            SetSampleDict();
        }

        /// <summary>
        /// Reset all values to null and lists to empty.
        /// </summary>
        public void Reset()
        {
            sequences.Clear();
            sequenceDict.Clear();
            sampleDict.Clear();
        }

        #endregion

        #region Private Methods


        // Todo fixme - creates far too many sequences
        private String GetSequenceDictKey(SAMAlignedSequence seq)
        {
            String seqStr = seq.QuerySequence.ToString();
            return Regex.Split(seqStr, "\r\n")[0];
        }


        // Todo fixme - results mostly look correct but don't match clstats
        private String GetSampleDictKey(SAMAlignedSequence seq)
        {
            foreach (SAMOptionalField field in seq.OptionalFields)
            {
                // I iterate through to find RG each time in case the optional fields
                // do not have a consistent format. Can I assume RG is always OptionalFields[0]
                // to avoid creating iterator and running this loop? todo
                if (field.Tag == "RG")
                {
                    return field.Value;
                }
            }
            return null;
        }

        /// <summary>
        /// Create a sequence dictionary from the current list of aligned sequences
        /// This is a dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value
        /// Each value in the dictionary is a shallow copy of an aligned sequence in sequences
        /// </summary>
        private void SetSequenceDict()
        {
            if (sequenceDict == null)
            {
                sequenceDict = new Dictionary<String, List<SAMAlignedSequence>>();
            }

            foreach (SAMAlignedSequence seq in sequences)
            {
                String sequenceStr = GetSequenceDictKey(seq);

                if (sequenceDict.ContainsKey(sequenceStr))
                {
                    List<SAMAlignedSequence> existingVal = sequenceDict[sequenceStr];
                    existingVal.Add(seq);
                    sequenceDict[sequenceStr] = existingVal;
                }
                else
                {
                    sequenceDict.Add(sequenceStr, new List<SAMAlignedSequence> { seq });
                }
            }
        }

        /// <summary>
        /// Create a sample dictionary from the current list of aligned sequences
        /// This is a dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// Each value in the dictionary is a shallow copy of an aligned sequence in sequences
        /// </summary>
        private void SetSampleDict()
        {
            if (sampleDict == null)
            {
                sampleDict = new Dictionary<String, List<SAMAlignedSequence>>();
            }

            foreach (SAMAlignedSequence seq in sequences)
            {
                String sampleName = GetSampleDictKey(seq);
                if (sampleName != null)
                {
                    if (sampleDict.ContainsKey(sampleName))
                    {
                        List<SAMAlignedSequence> existingVal = sampleDict[sampleName];
                        existingVal.Add(seq);
                        sampleDict[sampleName] = existingVal;
                    }
                    else
                    {
                        sampleDict.Add(sampleName, new List<SAMAlignedSequence> { seq });
                    }
                }
                else
                {
                    throw new ArgumentException("Sequence invalid: no sample identifier found");
                }
            }
        }

        #endregion

    }
}
