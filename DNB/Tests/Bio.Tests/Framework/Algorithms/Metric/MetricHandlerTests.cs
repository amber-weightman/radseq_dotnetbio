using System.Collections.Generic;
using System.Linq;
using Bio.Tests.Framework;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Bio.Algorithms.Metric;
using Bio.IO.SAM;
using System.Collections.ObjectModel;

namespace Bio.Tests.Algorithms.Alignment
{
    /// <summary>
    /// Tests for the MetricHandler classes.
    /// </summary>
    [TestClass]
    public class MetricHandlerTests
    {
        
        /// <summary>
        /// ...
        /// </summary>
        [TestMethod]
        [Priority(0)]
        [TestCategory("Priority0")]
        public void CMHNonDefaultConstructors()
        {
            string fileName = @"TestUtils\BAM\SeqAlignment.bam";
            Collection<SAMAlignedSequence> sequences = new Collection<SAMAlignedSequence>();
            sequences.Add(new SAMAlignedSequence());

            ClusterMetricHandler cmh = new ClusterMetricHandler(fileName);
            Assert.AreEqual(cmh.FileName, fileName);

            cmh = new ClusterMetricHandler(fileName, sequences);
            Assert.AreEqual(cmh.FileName, fileName);
            Assert.IsTrue(cmh.Sequences != null);
        }

        /// <summary>
        /// ...
        /// </summary>
        [TestMethod]
        [Priority(0)]
        [TestCategory("Priority0")]
        public void WriteToFile()
        {
            // file name is set
            // file type is set (sam/bam)
            // file type is valid
        }

        
        

    }
}
