using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Ploidulator
{
    /// <summary>
    /// Writes metric info to a csv file
    /// </summary>
    public class MetricFormatter : Bio.IO.IFormatter, IDisposable
    {
        #region Member variables

        /// <summary>
        /// Holds stream writer used for writing to file.
        /// </summary>
        private StreamWriter streamWriter = null;

        #endregion

        #region constructors

        /// <summary>
        /// Initializes a new instance of the MetricFormatter class.
        /// </summary>
        public MetricFormatter()
        {
            //throw new NotImplementedException();
        }

        /// <summary>
        /// Initializes a new instance of the MetricFormatter class with specified filename.
        /// </summary>
        /// <param name="fileName">MetricFormatter filename (including file extension).</param>
        public MetricFormatter(string fileName)
            : this()
        {
            this.Open(fileName);
        }
        #endregion

        #region Properties
        /// <summary>
        /// Gets the filename.
        /// </summary>
        public string FileName { get; private set; }

        /// <summary>
        /// Gets the name of this formatter.
        /// This is intended to give developers name of the formatter.
        /// </summary>
        public string Name
        {
            get
            {
                return Properties.Resources.METRIC_NAME;
            }
        }

        /// <summary>
        /// Gets the description of this formatter.
        /// This is intended to give developers some information 
        /// of the formatter class. This property returns a simple description of what this
        ///  class achieves.
        /// </summary>
        public string Description
        {
            get
            {
                return Properties.Resources.METRICFORMATTER_DESCRIPTION;
            }
        }

        /// <summary>
        /// Gets the file extension supported by this formatter.
        /// </summary>
        public string SupportedFileTypes
        {
            get
            {
                return Properties.Resources.METRIC_FILEEXTENSION;
            }
        }

        /// <summary>
        /// Gets or sets a value indicating whether the MetricFormatter will flush its buffer 
        /// to the underlying stream after every call to Write(ISequence).
        /// </summary>
        public bool AutoFlush { get; set; }

        #endregion

        #region Method
        /// <summary>
        /// Opens the specified file.
        /// </summary>
        /// <param name="fileName">Name of the file to open.</param>
        public void Open(string fileName)
        {
            if (this.streamWriter != null)
            {
                throw new InvalidOperationException(Properties.Resources.FILE_NOT_CLOSED);
            }

            this.FileName = fileName;
            this.streamWriter = new StreamWriter(this.FileName);
        }

        /// <summary>
        /// Opens the specified stream for writing sequences.
        /// </summary>
        /// <param name="outStream">StreamWriter to use.</param>
        public void Open(StreamWriter outStream)
        {
            if (this.streamWriter != null)
            {
                throw new InvalidOperationException(Properties.Resources.FILE_NOT_CLOSED);
            }

            this.FileName = null;
            this.streamWriter = outStream;
        }

        /// <summary>
        /// Writes an IMetric list as tab separated values to the file, one IMetric per line.
        /// </summary>
        /// <param name="metrics">IMetrics to write.</param>
        public void Write(IEnumerable<IMetric> metrics)
        {
            if (metrics == null)
            {
                throw new ArgumentNullException("metrics");
            }

            foreach (IMetric metric in metrics)
            {
                Write(metric);
            }

            this.streamWriter.Flush();
        }

        /// <summary>
        /// Writes the specified IMetric as tab separated values to a line of the file.
        /// </summary>
        /// <param name="metric">IMetric to write.</param>
        public void Write(IMetric metric)
        {
            if (metric == null)
            {
                throw new ArgumentNullException("metric");
            }

            if (this.streamWriter == null)
            {
                throw new InvalidOperationException(Properties.Resources.FILE_NOT_OPENED);
            }

            this.streamWriter.WriteLine(metric.ToString());

            if (this.AutoFlush)
            {
                this.streamWriter.Flush();
            }
        }

        /// <summary>
        /// Clears all buffer of underlying stream and any buffered data will be written to the file. 
        /// </summary>
        public void Flush()
        {
            if (this.streamWriter == null)
            {
                throw new InvalidOperationException(Properties.Resources.FILE_NOT_OPENED);
            }
            this.streamWriter.Flush();
        }

        /// <summary>
        /// Closes the current formatter and underlying stream.
        /// </summary>
        public void Close()
        {
            if (this.streamWriter == null)
            {
                throw new InvalidOperationException(Properties.Resources.FILE_NOT_OPENED);
            }

            this.Flush();
            this.streamWriter.Close();
            this.streamWriter.Dispose();
            this.streamWriter = null;
        }

        /// <summary>
        /// Disposes this formatter and underlying stream.
        /// </summary>
        public void Dispose()
        {
            if (this.streamWriter != null)
            {
                this.Close();
            }

            GC.SuppressFinalize(this);
        }

        #endregion
    }
}
