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

        /// <summary>
        /// Buffer used while writing to file.
        /// </summary>
        private byte[] buffer = null;

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
                return "Metric";
                //return Properties.Resource.METRIC_NAME;
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
                return "Writes an IMetric or List<IMetric> to a particular location, usually a file. The output is a tab separated text file with a .metr file extension";
                //return Properties.Resource.METRICFORMATTER_DESCRIPTION;
            }
        }

        /// <summary>
        /// Gets the file extension supported by this formatter.
        /// </summary>
        public string SupportedFileTypes
        {
            get
            {
                return ".metr";
                //return Properties.Resource.METRIC_FILEEXTENSION;
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
                throw new InvalidOperationException("File is already open.");
                //throw new InvalidOperationException(Bio.Properties.Resource.FileNotClosed);
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
                throw new InvalidOperationException("File is already open.");
                //throw new InvalidOperationException(Properties.Resource.FileNotClosed);
            }

            this.FileName = null;
            this.streamWriter = outStream;
        }


        // Writes a List<IMetric> as tab separated values to the file, one IMetric per line.

        // <param name="sequences">Sequences to write.</param>
        /*[Obsolete("Use the IEnumerable overload instead")]
        public void Write(ICollection<IMetric> metrics)
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
        }*/
        // todo fixme ok to remove?

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
                throw new InvalidOperationException("File is not opened. Please call Open method  to open the file.");
                //throw new InvalidOperationException(Properties.Resource.FileNotOpened);
            }

            string stringToWrite = metric.ToFileString();
            this.buffer = new byte[System.Text.ASCIIEncoding.Unicode.GetByteCount(stringToWrite)];
            this.streamWriter.WriteLine(metric.ToFileString());
            // todo fixme placeholder until we know what we are writing

            /*for (long index = 0; index < sequence.Count; index += maxLineSize)
            {
                for (bufferIndex = 0; bufferIndex < maxLineSize && index + bufferIndex < sequence.Count; bufferIndex++)
                {
                    this.buffer[bufferIndex] = sequence[index + bufferIndex];
                }

                string line = UTF8Encoding.UTF8.GetString(this.buffer, 0, bufferIndex);
                this.streamWriter.WriteLine(line);
            }*/

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
                throw new InvalidOperationException("File is not opened. Please call Open method  to open the file.");
                //throw new InvalidOperationException(Properties.Resource.FileNotOpened);
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
                throw new InvalidOperationException("File is not opened. Please call Open method  to open the file.");
                //throw new InvalidOperationException(Properties.Resource.FileNotOpened);
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



        // todo fixme reinstate this when i know what the format is
        /*public static string FormatString(ISequence sequence)
        {
            if (sequence == null)
            {
                throw new ArgumentNullException("sequence");
            }

            StringBuilder stringBuilder = new StringBuilder();

            stringBuilder.AppendLine(">" + sequence.ID);
            foreach (byte item in sequence)
            {
                stringBuilder.Append((char)item);
            }

            return stringBuilder.ToString();
        }*/

        #endregion
    }
}
