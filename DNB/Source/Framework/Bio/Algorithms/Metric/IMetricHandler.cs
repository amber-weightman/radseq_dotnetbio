using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio.Algorithms.Metric
{
    /// <summary>
    /// An IMetricHandler receives a set of aligned sequences or an individual sequence
    /// and performs operations on these sequence/s
    /// </summary>
    public interface IMetricHandler : IDisposable
    {
        /// <summary>
        /// Add a sequence.
        /// </summary>
        /// <param name="sequence">A sequence.</param>
        void Add(SAMAlignedSequence sequence);

        /// <summary>
        /// Add a list of sequences.
        /// </summary>
        /// <param name="sequences">A list of sequences.</param>
        void AddAll(List<SAMAlignedSequence> sequences);

        /// <summary>
        /// Process sequences (may process all, or a selection thereof).
        /// </summary>
        void ProcessSequences();

        /// <summary>
        /// Process any/all remaining sequences.
        /// </summary>
        void FlushSequences();

    }
}
