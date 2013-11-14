using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Ploidulator
{
    /// <summary>
    /// An IMetric calculates various metric values from a list of SAMAlignedSequences
    /// </summary>
    public interface IMetric
    {
        /// <summary>
        /// String of tab-separated values for writing to file by MetricFormatter.
        /// todo fixme we want the MetricFormater to decide how to write each line, but
        /// for now until I know what values are being written I will let the Metric handle this.
        /// </summary>
        string ToFileString();

        /// <summary>
        /// Calculate metric values from the given list of sequences.
        /// </summary>
        void Calculate(Collection<SAMAlignedSequence> clusterSequences);

        /// <summary>
        /// Reset all values to null and lists to empty.
        /// </summary>
        void Reset();

    }
}
