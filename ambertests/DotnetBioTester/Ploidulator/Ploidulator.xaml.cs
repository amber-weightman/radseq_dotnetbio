using Bio.Algorithms.Metric;
using Bio.IO.BAM;
using Bio.IO.SAM;
using Microsoft.Win32;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.DataVisualization;
using System.Windows.Controls.DataVisualization.Charting;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Animation;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Threading;

namespace Ploidulator
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private delegate System.Delegate QuickDelegate();
        private delegate System.Delegate IntDelegate(int a);
        private delegate System.Delegate HandlerDelegate(string a, string b, string c, string d, string d2,
            string e, string f, string f2,bool? g, bool? h, bool? i, bool? j, bool? k);
        private delegate System.Delegate StatsDelegate(int a, int b, int c, double d, double e, 
            double f, double g, double h, double i, double j, double k);

        private DispatcherTimer dispatcherTimer;
        private bool isProcessingFile = false;
        
        private ClusterMetricHandlerPloidulator handler = null;

        private DateTime startedTime;

        public MainWindow()
        {
            InitializeComponent();
        }

        #region button actions

        /// <summary>
        /// On click of GO button, begin ploidulation process
        /// </summary>
        private void Button_Click_Go(object sender, RoutedEventArgs e)
        {
            if (DataFileATextbox.Text != null)
            {
                // Launch ParseBAMMetric in a new thread
                HandlerDelegate handler = ParseBAMMetric;
                handler.BeginInvoke(DataFileATextbox.Text, ExpectedPloidyTextbox.Text, 
                    DirtCutoffTextbox.Text, AlignmentQualCutoffTextbox.Text,
                    ReadQualCutoffTextbox.Text, PopPercentTextbox.Text, GDirtCutoffTextbox.Text,
                    NumSamplesTextbox.Text, OutputToFile.IsChecked, 
                    MetricFileParent.IsChecked, MetricFileChild.IsChecked, OutputOverviewParent.IsChecked, 
                    OutputOverviewChild.IsChecked,
                    null, null);
            }
        }

        /// <summary>
        /// On click of X button, finish processing the current cluster and do not process any more
        /// </summary>
        private void Button_Click_Abort(object sender, RoutedEventArgs e)
        {
            // when the final handler thread finishes executing, it will hide this message
            AbortingMessage.Visibility = System.Windows.Visibility.Visible; 

            handler.Abort();
        }

        
        /// <summary>
        /// Open file selector window and fetch selected file path
        /// </summary>
        private void Button_Click_GetFile(object sender, RoutedEventArgs e)
        {
            OpenFileDialog dialog = new OpenFileDialog();
            dialog.DefaultExt = ".bam";
            dialog.Filter = "BAM files (*.bam)|*.bam";
            if ((bool)dialog.ShowDialog())
            {
                // Hide timer (may still be visible from previous ploidulation)
                TimerLabel.Visibility = System.Windows.Visibility.Hidden;

                // Display selected filename
                DataFileATextbox.Text = dialog.FileName;
            } 
        }

        #endregion

        #region timer

        /// <summary>
        /// Set/reset and display count-up timer
        /// </summary>
        private void StartTimer()
        {
            TimerLabel.Content = "";
            dispatcherTimer = new System.Windows.Threading.DispatcherTimer();
            dispatcherTimer.Tick += new EventHandler(dispatcherTimer_Tick);
            dispatcherTimer.Interval = new TimeSpan(0, 0, 1);
            dispatcherTimer.Start();
            startedTime = DateTime.Now;
            TimerLabel.Visibility = System.Windows.Visibility.Visible;
        }

        /// <summary>
        /// Pause/stop count-up timer
        /// </summary>
        private void StopTimer()
        {
            dispatcherTimer.Stop();
        }

        /// <summary>
        /// Update count-up timer label on dispatcher timer tick
        /// </summary>
        private void dispatcherTimer_Tick(object sender, EventArgs e)
        {
            TimeSpan time = DateTime.Now - startedTime;
            TimerLabel.Content = "Time elapsed: " + time.Hours.ToString() + ":" 
                + time.Minutes.ToString("00") + ":" + time.Seconds.ToString("00");
            // TODO is there a better C# way of doing this?
        }
        #endregion

        #region update gui

        /// <summary>
        /// Input file processing has begun. Update GUI to reflect this
        /// </summary>
        private System.Delegate UpdateGui_BeganParsing()
        {
            isProcessingFile = true;

            // Disable running any other new jobs
            ToggleSearchable(false);

            // Show stats panel
            StatsPanel.Visibility = System.Windows.Visibility.Visible;

            // Show progress bar and abort button
            ProgressStatusBar.Visibility = System.Windows.Visibility.Visible;
            AbortButton.Visibility = System.Windows.Visibility.Visible;

            // Show overview menus
            OverviewFileStats.IsExpanded = true;
            OverviewAllClustersStats.IsExpanded = true;
            OverviewGoodClustersStats.IsExpanded = false;

            StartTimer();

            return null;
        }

        /// <summary>
        /// Ploidulation for current input file has completed. Update GUI to reflect this
        /// </summary>
        /// <returns></returns>
        private System.Delegate UpdateGui_FinishedParsing()
        {
            // Enable running new job
            ToggleSearchable(true);

            // Hide progress bar
            ProgressStatusBar.Visibility = System.Windows.Visibility.Hidden;
            AbortButton.Visibility = System.Windows.Visibility.Hidden;

            // Hide any messages
            AbortingMessage.Visibility = System.Windows.Visibility.Hidden;

            // Reset to black
            SamplesDisplay.Foreground = Brushes.Black;
            SamplesDisplayLabel.Foreground = Brushes.Black;

            StopTimer();
            isProcessingFile = false;
            return null;
        }

        /// <summary>
        /// Enable or disable all search related fields
        /// </summary>
        private void ToggleSearchable(bool isSearchable)
        {
            FileExpander.IsExpanded = isSearchable;
            OutputFileExpander.IsExpanded = false;
            DataFileAButton.IsEnabled = isSearchable;
            DataFileATextbox.IsEnabled = isSearchable;
            ExpectedPloidyTextbox.IsEnabled = isSearchable;
            AlignmentQualCutoffTextbox.IsEnabled = isSearchable;
            DirtCutoffTextbox.IsEnabled = isSearchable;
            NumSamplesTextbox.IsEnabled = isSearchable;
            ReadQualCutoffTextbox.IsEnabled = isSearchable;
            PopPercentTextbox.IsEnabled = isSearchable;
            LaunchButton.IsEnabled = isSearchable;
        }

        /// <summary>
        /// Update stats panels with data from handler
        /// </summary>
        private System.Delegate UpdateStatsPanel(int numGoodClusters, int numClustersParsed, int maxSampleCount,
            double maxMapq, double maxReadq, double avgDirt, double avgMapq, double avgReadq, 
            double avgDirtGood, double avgMapqGood, double avgReadqGood)
        {
            if (numClustersParsed != 0)
            {
                UpdateProgress(numClustersParsed);

                SamplesDisplay.Content = maxSampleCount.ToString();
                MAPQdisplay.Content = maxMapq.ToString();
                READQdisplay.Content = maxReadq.ToString();

                DirtDisplay.Content = avgDirt.ToString();
                AvgMAPQDisplay.Content = avgMapq.ToString();
                AvgREADQDisplay.Content = avgReadq.ToString();

                DirtDisplayGood.Content = avgDirtGood.ToString();
                AvgMAPQDisplayGood.Content = avgMapqGood.ToString();
                AvgREADQDisplayGood.Content = avgReadqGood.ToString();

                double numGood = Math.Round(numGoodClusters / (double)numClustersParsed * 100, 2);
                double numBad = 100 - numGood;
                GoodBadProgressBar.Value = numGood;
                GoodCountLabel.Content = numGood + "% good";
                BadCountLabel.Content = numBad + "% bad";

                string plural = numGoodClusters == 1 ? "" : "s";
                string plural2 = (numClustersParsed - numGoodClusters) == 1 ? "" : "s";
                GoodCountLabel.ToolTip = numGood + "% of clusters are good (" + numGoodClusters + " cluster" + plural + ")";
                BadCountLabel.ToolTip = numBad + "% of clusters are bad (" + (numClustersParsed - numGoodClusters) + " cluster" + plural2 + ")";

                int a = Convert.ToInt32(maxSampleCount);
                int b = Convert.ToInt32(NumSamplesTextbox.Text);
                if ((int)a < (int)b && (int)a != 0)
                {
                    SamplesDisplay.Foreground = Brushes.Red;
                    SamplesDisplayLabel.Foreground = Brushes.Red;
                }
                else
                {
                    SamplesDisplay.Foreground = Brushes.Black;
                    SamplesDisplayLabel.Foreground = Brushes.Black;
                }

                if (numGoodClusters > 0)
                {
                    OverviewGoodClustersStats.IsExpanded = true;
                }

            }
            return null;
        }

        /// <summary>
        /// Set initial progress bar values (0 to max) based on number of clusters in input file
        /// </summary>
        /// <param name="max"></param>
        private System.Delegate SetProgressBar(int max)
        {
            ProgressBar.IsEnabled = true;
            ProgressBar.Minimum = 0;
            ProgressBar.Maximum = max;

            if(max != 0)
            {    
                ProgressBar.IsIndeterminate = false;
                ProgressBar.ToolTip = "Found " + max.ToString() + " clusters to analyse ...";
            }
            else
            {
                ProgressBar.IsIndeterminate = true;
                ProgressBar.ToolTip = "Processing clusters...";
            }
            ProgressBar.Visibility = System.Windows.Visibility.Visible;
            AbortButton.Visibility = System.Windows.Visibility.Visible;
            LoadingBarLabel.Visibility = System.Windows.Visibility.Visible;
                
            return null;
        }

        /// <summary>
        /// Update the progress bar value
        /// </summary>
        private System.Delegate UpdateProgress(int index)
        {
            if (index >= ProgressBar.Minimum && index <= ProgressBar.Maximum)
            {
                ProgressBar.Value = index;
                double percent = Math.Round((index / (double)ProgressBar.Maximum * 100), 2);
                ProgressBar.ToolTip = "Analysed " + index.ToString() + " out of " + ProgressBar.Maximum.ToString() + " clusters (" + percent.ToString() + "%)";
                LoadingBarLabel.Content = ProgressBar.ToolTip;
            }
            else if (ProgressBar.Maximum == 0)
            {
                ProgressBar.ToolTip = "Analysed " + index.ToString() + " clusters";
                LoadingBarLabel.Content = ProgressBar.ToolTip;
            }
            return null;
        }

        #endregion

        #region process sequences

        /// <summary>
        /// Initialise handler and handler settings based on user's input
        /// </summary>
        private void InitHandler(string filename, string ploidy, string dirtCutoff, string alignQualCutoff,
            string readQualCutoff, string popPercent, string gDirtCutoff, string numSamples, bool? outputToFile, bool? metricFileParent,
            bool? metricFileChild, bool? outputOverviewParent, bool? outputOverviewChild)
        {
            // Check user input
            string newName = filename.Split(new char[] { '.' })[0];

            // Set up the handler
            handler = new ClusterMetricHandlerPloidulator(newName + "_pl", Convert.ToInt32(ploidy), Convert.ToDouble(dirtCutoff),
                Convert.ToDouble(alignQualCutoff), Convert.ToDouble(readQualCutoff), Convert.ToDouble(popPercent), Convert.ToInt32(numSamples), outputToFile);

            handler.GDirtCutoff = Convert.ToDouble(gDirtCutoff);
            handler.WriteToFilteredBam = outputToFile == true;
            handler.WriteClusterMetricOriginal = metricFileParent == true;
            handler.WriteClusterMetricFiltered = metricFileChild == true;
            handler.WriteOverviewMetricOriginal = outputOverviewParent == true;
            handler.WriteOverviewMetricFiltered = outputOverviewChild == true;
        }

        /// <summary>
        /// Read input file header and set progress bar based on number of sequences in input file
        /// </summary>
        private void ReadHeader(string filename, BAMParser parser, int numClustersInInputFile)
        {
            using (Stream readStream = new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                handler.InputHeader = parser.GetHeader(readStream);
                numClustersInInputFile = handler.InputHeader.ReferenceSequences.Count;
                Console.WriteLine("header has clust count " + numClustersInInputFile);
                /*foreach(ReferenceSequenceInfo cluster in header.ReferenceSequences)
                {
                    Console.Write(cluster.Name + " - ");
                }*/
                //header.RecordFields["RG"];
            }
            Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
                new IntDelegate(SetProgressBar), numClustersInInputFile);
        }

        /// <summary>
        /// Create handler and initialise handler settings, then process sequences
        /// </summary>
        private System.Delegate ParseBAMMetric(string filename, string ploidy, string dirtCutoff, string alignQualCutoff, 
            string readQualCutoff, string popPercent, string gDirtCutoff, string numSamples, bool? outputToFile, 
            bool? metricFileParent, bool? metricFileChild, bool? outputOverviewParent, bool? outputOverviewChild)
        {
            BAMParser parser = new BAMParser();
            int numClustersInInputFile = 0;
            
            // Initialise the metric handler
            InitHandler(filename, ploidy, dirtCutoff, alignQualCutoff, readQualCutoff, popPercent, gDirtCutoff, numSamples, 
                outputToFile, metricFileParent, metricFileChild, outputOverviewParent, outputOverviewChild);

            ReadHeader(filename, parser, numClustersInInputFile);

            // Begin parsing sequences
            Dispatcher.BeginInvoke( System.Windows.Threading.DispatcherPriority.Normal,
                new QuickDelegate(UpdateGui_BeganParsing)); // updates GUI to indicate that parsing has begun
            ProcessSequences(filename, parser);
                
            // Finished sequences for this input file
            Dispatcher.BeginInvoke( System.Windows.Threading.DispatcherPriority.Normal,
                new QuickDelegate(UpdateGui_FinishedParsing));

            // Close handler and return
            Console.WriteLine("FINISHED");
            handler.Dispose();
            parser.Dispose();
            return null;
        }

        /// <summary>
        /// For each sequence read by the parser, pass it to the metric handler. Also update GUI summary panels
        /// </summary>
        private void ProcessSequences(string filename, BAMParser parser)
        {
            double incrementOne = 1; // the progress bar will update every [incrementOne] many clusters
            double incrementTwo = 1; // the progress bar will update every [incrementTwo] many clusters
            int increaseIncrementAfter = 100; // after this many clusters, increase the increment at which we update the gui
            double increment;       // the current increment at which progress bar will update
            int updateDisplayForClusterIndex = -1, clusterCount = -1; // whether we have already updated the gui for this cluster
            foreach (SAMAlignedSequence se in parser.ParseSequence(filename))
            {
                if (handler.Add(se))
                {
                    // optionalFields["RG"]
                    clusterCount = handler.ClusterCount;
                    increment = (clusterCount < increaseIncrementAfter) ? incrementOne : incrementTwo;
                    if ((clusterCount == 1 || (double)clusterCount % increment == 0) && clusterCount != updateDisplayForClusterIndex)
                    {
                        updateDisplayForClusterIndex = clusterCount;

                        Dispatcher.BeginInvoke(
                            System.Windows.Threading.DispatcherPriority.Normal,
                            new StatsDelegate(UpdateStatsPanel), handler.GoodCount, clusterCount,
                                handler.MaxSampleCount, handler.MaxMapQ, handler.MaxReadQ,
                                handler.AverageDirt, handler.AverageMapQ, handler.AverageReadQ,
                                handler.AverageDirtGood, handler.AverageMapQGood, handler.AverageReadQGood);
                    }
                }
                else
                {
                    Console.WriteLine("Handler has finished processing sequences");
                    break;
                }
            }
            // Tell the handler that there are no more sequences to receive
            handler.SetComplete(); 
        }
      
        #endregion

        #region window and controls

        /// <summary>
        /// If user attempts to close GUI window while data is still being written to file,
        /// give warning
        /// </summary>
        private void Window_Closing(object sender, System.ComponentModel.CancelEventArgs e)
        {
            if(isProcessingFile)
            {
                MessageBoxResult result = MessageBox.Show("Data is currently being output to file. \nIf you close this window the output file will be incomplete.\n\nAre you sure you want to close this window?",
                "Warning", MessageBoxButton.YesNo, MessageBoxImage.Warning);
                if (result == MessageBoxResult.No)
                {
                    e.Cancel = true;
                }
            }
        }

        /// <summary>
        /// Prevent non-numeric input in a textbox
        /// </summary>
        private void NumericOnly_TextChanged(object sender, TextChangedEventArgs e)
        {
            // todo fixme - multiple ways of doing this, this is probably not the most efficient
            TextBox t = sender as TextBox;
            string s = "";
            int dec = 0;
            foreach (Char c in t.Text.ToCharArray())
            {
                s += Char.IsDigit(c) || (c == '.' && dec == 0)? c.ToString() : "";
                if (c == '.')
                {
                    dec++;
                }
            }
            t.Text = s;
        }

        #endregion

        #region visualisation

        /// <summary>
        /// Draws a piechart. Not currently in use
        /// </summary>
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
            wpMain.Children.Add(c);
            return null;
        }

        #endregion

    }
}
