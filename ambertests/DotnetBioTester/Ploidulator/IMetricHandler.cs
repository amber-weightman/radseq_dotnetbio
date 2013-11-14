using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Ploidulator
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
        /// <returns>True if successfully added, false otherwise</returns>
        bool Add(SAMAlignedSequence sequence);

        /// <summary>
        /// Add a set of sequences.
        /// THIS IS NOT ACTUALLY CALLED IN BAMPARSER, BUT SEEMS A LOGICAL INCLUSION - AW
        /// </summary>
        /// <param name="sequences">A set of sequences.</param>
        /// <returns>True if successfully added, false otherwise</returns>
        bool AddRange(IEnumerable<SAMAlignedSequence> sequences);

        /// <summary>
        /// Process sequences (may process all, or a selection thereof).
        /// ALSO NOT CALLED BY BAMPARSER. IS THIS A LOGICAL INCLUSION FOR A MORE FUNCTIONAL
        /// INTERFACE, OR SHOULD IT BE REMOVED? - AW
        /// </summary>
        void ProcessSequences();

        /// <summary>
        /// Process any/all remaining sequences.
        /// </summary>
        void SetComplete();

    }
}
