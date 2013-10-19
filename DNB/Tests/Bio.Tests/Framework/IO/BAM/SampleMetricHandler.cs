using Bio.Algorithms.Metric;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;

namespace Bio.Tests.Framework.IO.BAM
{
    /// <summary>
    /// Metric handler shell, for testing.
    /// No assumptions should be made about how sequences are handled. This sample class
    /// stores all sequences added in Sequences. ProcessSequences removes a single sequence, 
    /// and FlushSequences removes all sequences
    /// </summary>
    class SampleMetricHandler : IMetricHandler
    {

        public SampleMetricHandler()
        {
            sequences = new List<SAMAlignedSequence>();
        }

        /// <summary>
        /// Sequences stored, not yet processed.
        /// </summary>
        private List<SAMAlignedSequence> sequences = null;

        /// <summary>
        /// Sequences stored, not yet processed.
        /// </summary>
        public List<SAMAlignedSequence> Sequences
        {
            get { return sequences; }
        }

        /// <summary>
        /// Add a sequence.
        /// </summary>
        /// <param name="sequence">Sequence to add.</param>
        public void Add(SAMAlignedSequence sequence)
        {
            sequences.Add(sequence);
        }

        /// <summary>
        /// Add a list of sequences.
        /// </summary>
        /// <param name="seqs">Sequences to add.</param>
        public void AddRange(IEnumerable<SAMAlignedSequence> seqs)
        {
            foreach(SAMAlignedSequence seq in seqs)
            {
                sequences.Add(seq);
            }
        }

        /// <summary>
        /// Remove one sequence from the stored list.
        /// </summary>
        public void ProcessSequences()
        {
            sequences.RemoveAt(0);
        }

        /// <summary>
        /// Remove all sequences from the stored list.
        /// </summary>
        public void SetComplete()
        {
            sequences = new List<SAMAlignedSequence>();
        }

        /// <summary>
        /// Dispose method.
        /// </summary>
        public void Dispose()
        {
            //throw new NotImplementedException();
        }
    }
}
