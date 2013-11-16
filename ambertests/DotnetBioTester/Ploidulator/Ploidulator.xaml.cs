using Accord.Statistics.Distributions.Univariate;
using Bio.IO.BAM;
using Bio.IO.SAM;
using Microsoft.Win32;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.DataVisualization.Charting;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Threading;

namespace Ploidulator
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml, being the Ploidulator user interface
    /// </summary>
    public partial class MainWindow : Window
    {
        private static CultureInfo ci = new CultureInfo("en-AU");

        #region private fields

        #region delegates
        private delegate System.Delegate QuickDelegate();
        private delegate System.Delegate IntDelegate(int a);
        private delegate System.Delegate SingleIntChartDataDelegate(int c, int d);
        private delegate System.Delegate ChartTupleDelegate(IEnumerable a, IEnumerable b);
        private delegate System.Delegate ChartDelegate(Chart a);
        private delegate System.Delegate IntChartDataDelegate(KeyValuePair<int, int>[] a, KeyValuePair<int, int>[] b);
        private delegate System.Delegate HandlerDelegate(string a, string b, string c, string d, string d1, string d2,
            string e, string f, bool? f1, string f2, bool? g, bool? h, bool? i, bool? j, bool? k, bool? l, bool? m);
        private delegate System.Delegate StatsDelegate(int a, int b, int c, double d, double e, 
            double f, double g, double h, double i, double j, double k);
        #endregion

        #region table data

        /// <summary>
        /// Data series A for columnchart (read distribution - all distinct reads)
        /// </summary>
        private KeyValuePair<int, double>[] columnSeriesA = null;

        /// <summary>
        /// Data series B for columnchart (read distribution - good distinct reads)
        /// </summary>
        private KeyValuePair<int, double>[] columnSeriesB = null;

        #endregion

        #region timers
        private DispatcherTimer TimeElapsedTimer;   // Timer records time elapsed while a file is being parsed
        private DispatcherTimer ChartUpdateTimer;   // Charts are updated at regular intervals on a timer
        private DateTime startedTime;
        #endregion

        #region flags
        private bool isProcessingFile = false;      // Whether an input file is currently being parsed
        #endregion

        private ClusterMetricHandler metricHandler = null; // Handler object for parsing the input file

        #endregion

        #region public methods

        public MainWindow()
        {
            InitializeComponent();
        }

        #endregion

        #region private methods

        #region buttons

        /// <summary>
        /// On click of GO button, begin ploidulation process
        /// </summary>
        private void Button_Click_Go(object sender, RoutedEventArgs e)
        {
            // Input file and number of samples must not be null
            if (!String.IsNullOrEmpty(DataFileATextbox.Text) &&
                NumSamplesTextbox.Text != null && NumSamplesTextbox.Text != "0"
                && File.Exists(DataFileATextbox.Text))
            {
                // Launch ParseBAMMetric in a new thread
                HandlerDelegate handler = ParseBAMMetric;
                handler.BeginInvoke(DataFileATextbox.Text, ExpectedPloidyTextbox.Text,
                    DirtCutoffTextbox.Text, PloidyDisagreementTextbox.Text, AlignmentQualCutoffTextbox.Text,
                    ReadQualCutoffTextbox.Text, PopPercentTextbox.Text, MaxHaplotypesTextbox.Text, OnlyHaplotypeGood.IsChecked,
                    NumSamplesTextbox.Text, OutputToFile.IsChecked, 
                    MetricFileParent.IsChecked, MetricFileChild.IsChecked, OutputOverviewParent.IsChecked, 
                    OutputOverviewChild.IsChecked, GenotypesToFile.IsChecked, HaplotypesToFile.IsChecked, 
                    null, null);
            }
            else
            {
                MessageBox.Show(Properties.Resources.MANDATORY_FIELDS_WARNING,
                Properties.Resources.MANDATORY_FIELDS_WARNING_HEADER, MessageBoxButton.OK, MessageBoxImage.Asterisk);
            }
        }

        /// <summary>
        /// On click of X button, finish processing the current cluster and do not process any more
        /// </summary>
        private void Button_Click_Abort(object sender, RoutedEventArgs e)
        {
            if(metricHandler != null)
            {
                // when the final handler thread finishes executing, it will hide this message
                AbortingMessage.Visibility = System.Windows.Visibility.Visible; 
                metricHandler.Abort();
            }
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
            CreateTimer(ref TimeElapsedTimer, new EventHandler(dispatcherTimer_Tick), 0, 0, 1);
            startedTime = DateTime.Now;
            TimerLabel.Visibility = System.Windows.Visibility.Visible;
        }

        /// <summary>
        /// Set/reset and chart updating timer
        /// (charts are updated on a timer to avoid rendering too frequently if sequences are parsed quickly)
        /// </summary>
        private void StartChartTimer()
        {
            CreateTimer(ref ChartUpdateTimer, new EventHandler(chartTimer_Tick), 0, 0, 10); // update every 10 seconds
        }

        /// <summary>
        /// Create a timer
        /// </summary>
        private static void CreateTimer(ref DispatcherTimer timer, EventHandler eventHandler, int h, int m, int s)
        {
            timer = new System.Windows.Threading.DispatcherTimer();
            timer.Tick += eventHandler;
            timer.Interval = new TimeSpan(h, m, s);
            timer.Start();
        }

        /// <summary>
        /// Pause/stop all timers
        /// </summary>
        private void StopTimers()
        {
            TimeElapsedTimer.Stop();
            ChartUpdateTimer.Stop();
        }

        /// <summary>
        /// Update count-up timer label on dispatcher timer tick
        /// </summary>
        private void dispatcherTimer_Tick(object sender, EventArgs e)
        {
            TimeSpan time = DateTime.Now - startedTime;
            TimerLabel.Content = "Time elapsed: " + time.Hours.ToString(ci) + ":" 
                + time.Minutes.ToString("00",ci) + ":" + time.Seconds.ToString("00",ci);
        }

        /// <summary>
        /// Update displayed chart/s on timer tick
        /// </summary>
        private void chartTimer_Tick(object sender, EventArgs e)
        {
           UpdateDisplayedChart();
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
            AbortMenuItem.IsEnabled = true;

            // Show overview menus
            OverviewFileStats.IsExpanded = true;
            OverviewAllClustersStats.IsExpanded = true;
            OverviewGoodClustersStats.IsExpanded = false;

            Tab1.IsEnabled = true;
            Tab2.IsEnabled = true;
            Tab3.IsEnabled = true;

            StartTimer();
            StartChartTimer();

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
            ProgressBar.Value = 0;
            LoadingBarLabel.Content = "";
            AbortButton.Visibility = System.Windows.Visibility.Hidden;
            AbortMenuItem.IsEnabled = false;

            // Hide any messages
            AbortingMessage.Visibility = System.Windows.Visibility.Hidden;

            // Reset to black
            SamplesDisplay.Foreground = Brushes.Black;
            SamplesDisplayLabel.Foreground = Brushes.Black;

            StopTimers();
            isProcessingFile = false;
            return null;
        }

        /// <summary>
        /// Enable or disable all search related fields
        /// </summary>
        private void ToggleSearchable(bool isSearchable)
        {
            DataFileAButton.IsEnabled = isSearchable;
            DataFileATextbox.IsEnabled = isSearchable;
            ExpectedPloidyTextbox.IsEnabled = isSearchable;
            AlignmentQualCutoffTextbox.IsEnabled = isSearchable;
            DirtCutoffTextbox.IsEnabled = isSearchable;
            NumSamplesTextbox.IsEnabled = isSearchable;
            ReadQualCutoffTextbox.IsEnabled = isSearchable;
            PopPercentTextbox.IsEnabled = isSearchable;
            MaxHaplotypesTextbox.IsEnabled = isSearchable;
            OnlyHaplotypeGood.IsEnabled = isSearchable;
            LaunchButton.IsEnabled = isSearchable;
            MetricFileParent.IsEnabled = isSearchable;
            MetricFileChild.IsEnabled = isSearchable;
            OutputOverviewParent.IsEnabled = isSearchable;
            OutputOverviewChild.IsEnabled = isSearchable;
            OutputToFile.IsEnabled = isSearchable;
            PloidyDisagreementTextbox.IsEnabled = isSearchable;
            GenotypesToFile.IsEnabled = isSearchable;
            HaplotypesToFile.IsEnabled = isSearchable;
        }

        
        /// <summary>
        /// Update stats panels with data from handler
        /// </summary>
        private Delegate UpdateStatsPanel(int numGoodClusters, int numClustersParsed, int maxSampleCount,
            double maxMapq, double maxReadq, double avgDirt, double avgMapq, double avgReadq, 
            double avgDirtGood, double avgMapqGood, double avgReadqGood)
        {
            if (numClustersParsed != 0)
            {
                UpdateProgress(numClustersParsed);

                SamplesDisplay.Content = maxSampleCount.ToString(ci);
                MAPQdisplay.Content = maxMapq.ToString(ci);
                READQdisplay.Content = maxReadq.ToString(ci);

                DirtDisplay.Content = avgDirt.ToString(ci);
                AvgMAPQDisplay.Content = avgMapq.ToString(ci);
                AvgREADQDisplay.Content = avgReadq.ToString(ci);

                // If this is the first time we have had good clusters, expand the good clusters menu
                if (numGoodClusters > 0 && !String.IsNullOrEmpty(DirtDisplayGood.Content.ToString()))
                {
                    OverviewGoodClustersStats.IsExpanded = true;
                }

                DirtDisplayGood.Content = avgDirtGood.ToString(ci);
                AvgMAPQDisplayGood.Content = avgMapqGood.ToString(ci);
                AvgREADQDisplayGood.Content = avgReadqGood.ToString(ci);

                double numGood = Math.Round(numGoodClusters / (double)numClustersParsed * 100, 2);
                double numBad = 100 - numGood;
                GoodBadProgressBar.Value = numGood;
                GoodCountLabel.Content = numGood + "% good";
                BadCountLabel.Content = numBad + "% bad";

                string plural = numGoodClusters == 1 ? "" : "s";
                string plural2 = (numClustersParsed - numGoodClusters) == 1 ? "" : "s";
                GoodCountLabel.ToolTip = numGood + "% of clusters are good (" + numGoodClusters + " cluster" + plural + ")";
                BadCountLabel.ToolTip = numBad + "% of clusters are bad (" + (numClustersParsed - numGoodClusters) + " cluster" + plural2 + ")";

                int a = Convert.ToInt32(maxSampleCount, ci);
                int b = Convert.ToInt32(NumSamplesTextbox.Text, ci);
                if ((int)a != (int)b && (int)a != 0)
                {
                    SamplesDisplay.Foreground = Brushes.Red;
                    SamplesDisplayLabel.Foreground = Brushes.Red;
                }
                else
                {
                    SamplesDisplay.Foreground = Brushes.Black;
                    SamplesDisplayLabel.Foreground = Brushes.Black;
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
                ProgressBar.ToolTip = "Found " + max.ToString(ci) + " clusters to analyse ...";
            }
            else
            {
                ProgressBar.IsIndeterminate = true;
                ProgressBar.ToolTip = "Processing clusters...";
            }
            ProgressBar.Visibility = System.Windows.Visibility.Visible;
            AbortButton.Visibility = System.Windows.Visibility.Visible;
            LoadingBarLabel.Visibility = System.Windows.Visibility.Visible;
            AbortMenuItem.IsEnabled = true;
                
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
                ProgressBar.ToolTip = "Analysed " + index.ToString(ci) + " out of " + ProgressBar.Maximum.ToString(ci) + " clusters (" + percent.ToString(ci) + "%)";
                LoadingBarLabel.Content = ProgressBar.ToolTip;
            }
            else if (ProgressBar.Maximum == 0)
            {
                ProgressBar.ToolTip = "Analysed " + index.ToString(ci) + " clusters";
                LoadingBarLabel.Content = ProgressBar.ToolTip;
            }
            return null;
        }


        #endregion

        #region visualisation

        /// <summary>
        /// Update the ItemsSource for both pie charts
        /// </summary>
        private Delegate UpdatePiechart(IEnumerable itemsSourceA, IEnumerable itemsSourceB)
        {
            ((PieSeries)PloidyPie.Series[0]).ItemsSource = itemsSourceA;
            if (itemsSourceB != null)
            {
                ((PieSeries)PloidyPieGood.Series[0]).ItemsSource = itemsSourceB;
                PloidyPieGood.Visibility = System.Windows.Visibility.Visible;
            }
            return null;
        }

        /// <summary>
        /// Returns a tuple of pie chart ItemsSources (enables partial rendering of piechart in background thread)
        /// </summary>
        private static Tuple<IEnumerable, IEnumerable> PreRenderPiechart(KeyValuePair<string, double>[] dataSeriesA, 
            KeyValuePair<string, double>[] dataSeriesB)
        {
            Chart temp = new Chart();
            DrawPieChart(ref temp, dataSeriesA);
            if (dataSeriesB != null)
            {
                Chart tempB = new Chart();
                DrawPieChart(ref tempB, dataSeriesB);
                return new Tuple<IEnumerable, IEnumerable>(((PieSeries)temp.Series[0]).ItemsSource, 
                    ((PieSeries)tempB.Series[0]).ItemsSource);
            }
            return new Tuple<IEnumerable, IEnumerable>(((PieSeries)temp.Series[0]).ItemsSource, null);
        }

        /// <summary>
        /// Update the ItemsSource for the column chart
        /// </summary>
        private Delegate UpdateColumnchart(int min, int max)
        {
            //((ColumnSeries)ReadCountChart.Series[0]).ItemsSource = null;
            ((ColumnSeries)ReadCountChart.Series[0]).ItemsSource = columnSeriesA;
            if (columnSeriesB != null && columnSeriesB.Length > 0)
            {
                //((ColumnSeries)ReadCountChart.Series[1]).ItemsSource = null;
                ((ColumnSeries)ReadCountChart.Series[1]).ItemsSource = columnSeriesB;
            }
            //DrawPoisson(chart, dataSeries);

            SetColumnchartSliders(0, max);
            return null;
        }

        /// <summary>
        /// Set the min and max values for the columnchart zoom in/out sliders
        /// </summary>
        private void SetColumnchartSliders(int min, int max)
        {
            // Range covered by both sliders must be between min and max
            SliderA.Minimum = min;
            SliderB.Maximum = max;

            // If the second slider value has not yet been changed by the user, set the selected value to max
            // (zooms slder out)
            if (SliderB.SelectionEnd == 0) 
            {
                SliderA.Maximum = max-1;
                SliderB.Minimum = 1;
                SliderB.Value = max;
            }
        }

        /// <summary>
        /// Update the ItemsSource for the linechart
        /// </summary>
        private Delegate UpdateLinechart(Chart chart)
        {
            ((LineSeries)IndivChart.Series[0]).ItemsSource = ((LineSeries)chart.Series[0]).ItemsSource;
            if (chart.Series.Count > 1)
            {
                ((LineSeries)IndivChart.Series[1]).ItemsSource = ((LineSeries)chart.Series[1]).ItemsSource;
            } 
            return null;
        }

        /// <summary>
        /// Update whichever chart is currently displayed in the active tab
        /// </summary>
        private void UpdateDisplayedChart()
        {
            if (Tab1.IsSelected)
            {
                UpdatePiechartDisplay(null, null);
            }
            else if (Tab2.IsSelected)
            {
                UpdateColumnchartDisplay(null, null);
            }
            else if (Tab3.IsSelected)
            {
                UpdateLinechartDisplay(null, null);
            }
        }

        
        /// <summary>
        /// Update the piechart
        /// </summary>
        private void UpdatePiechartDisplay(object sender, RoutedEventArgs e)
        {
            int ploidyLevel = Convert.ToInt32(ExpectedPloidyTextbox.Text, ci);

            KeyValuePair<string, double>[] dataSeriesA = GetPieData(metricHandler.ClusterSequenceFrequenciesOverview, ploidyLevel);
            KeyValuePair<string, double>[] dataSeriesB = GetPieData(metricHandler.ClusterSequenceFrequenciesOverviewGood, ploidyLevel);

            Tuple<IEnumerable, IEnumerable> charts = PreRenderPiechart(dataSeriesA, dataSeriesB);

            Dispatcher.BeginInvoke( System.Windows.Threading.DispatcherPriority.Normal,
                new ChartTupleDelegate(UpdatePiechart), charts.Item1, charts.Item2);
        }

        /// <summary>
        /// Update the columnchart
        /// </summary>
        private void UpdateColumnchartDisplay(object sender, RoutedEventArgs e)
        {
            double rounding = Convert.ToDouble(ReadCountRounding.Text, ci);
            //KeyValuePair<int, double>[] data = GetColumnData(metricHandler.GraphDataDistinctReads, metricHandler.NumberClustersParsed, rounding);
            //KeyValuePair<int, double>[] datab = GetColumnData(metricHandler.GraphDataDistinctReadsGood, metricHandler.GoodCount, rounding);
            columnSeriesA = GetColumnData(metricHandler.GraphDataDistinctReads, metricHandler.NumberClustersParsed, rounding);
            columnSeriesB = GetColumnData(metricHandler.GraphDataDistinctReadsGood, metricHandler.GoodCount, rounding);

            int min = metricHandler.MinDistinctReadCount;
            int max = metricHandler.MaxDistinctReadCount;

            Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
                new SingleIntChartDataDelegate(UpdateColumnchart), min, max);
        }

        /// <summary>
        /// retrieve data for the column chart, rounded by [rounding] and sorted
        /// </summary>
        private static KeyValuePair<int, double>[] GetColumnData(Dictionary<int, int> data, double count, double rounding)
        {
            Dictionary<int, double> roundedData = new Dictionary<int, double>();
            foreach (KeyValuePair<int, int> entry in data)
            {
                double val = Math.Round((double)(entry.Value / (double)count * 100), 2); // Percentage of clusters, to nearest two decimals
                int key = (int)Math.Round(entry.Key / rounding, 0) * (int)rounding; ; // Number of reads (rounded to nearest n)
                if (key >= rounding)
                {
                    if (roundedData.ContainsKey(key))
                    {
                        roundedData[key] += val;
                    }
                    else
                    {
                        roundedData[key] = val;
                    }
                }
            }
            var sorted = (from datum in roundedData orderby datum.Key ascending select datum)
                    .ToDictionary(pair => pair.Key, pair => pair.Value);
            return sorted.ToArray();
        }

        /// <summary>
        /// Update the linechart
        /// </summary>
        private void UpdateLinechartDisplay(object sender, RoutedEventArgs e)
        {
            Chart chart = PreDrawLineChart(metricHandler.GraphDataIndividualsCounts, metricHandler.GraphDataIndividualsCountsGood);
            Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal, new ChartDelegate(UpdateLinechart), chart);
        }

        /// <summary>
        /// Draw a Poisson curve over the columnchart
        /// </summary>
        [System.Diagnostics.CodeAnalysis.SuppressMessage("Microsoft.Performance", "CA1811:AvoidUncalledPrivateCode")]
        private void DrawPoisson(Chart chart, KeyValuePair<int, double>[] dataSeries)
        {
            // Calculate Poisson mean
            int max = dataSeries[dataSeries.Length - 1].Key;
            double mean = metricHandler.ReadCountDistinctTotal / max;

            PoissonDistribution dist = new PoissonDistribution(lambda: mean);
            List<KeyValuePair<int, double>> dataSeries2 = new List<KeyValuePair<int, double>>();
            int i = 1;
            double pdf = 0;
            while (true)
            {
                pdf = dist.ProbabilityMassFunction(k: i);
                if (pdf == 0 || Double.IsNaN(pdf)) { break; }
                double val = pdf * 100;
                dataSeries2.Add(new KeyValuePair<int, double>(i, val));
                i += 2;
            }
            ((LineSeries)chart.Series[2]).ItemsSource = dataSeries2;
        }

        /// <summary>
        /// Enables the linechart to be partially rendered in the background thread
        /// </summary>
        private static Chart PreDrawLineChart(Dictionary<int, int> dataSeriesA, Dictionary<int, int> dataSeriesB)
        {
            Chart chart = new Chart();
            if (chart.Series.Count == 0)
            {
                chart.Series.Add(new LineSeries());
            }
            ((LineSeries)chart.Series[0]).ItemsSource = dataSeriesA;
            if (dataSeriesB != null && dataSeriesB.Count > 0)
            {
                if (chart.Series.Count == 1)
                {
                    chart.Series.Add(new LineSeries());
                }
                ((LineSeries)chart.Series[1]).ItemsSource = dataSeriesB;
            }
            return chart;
        }

        /// <summary>
        /// Draw a piechart
        /// </summary>
        private static void DrawPieChart(ref Chart chart, KeyValuePair<string, double>[] dataSeries)
        {
            if(dataSeries == null || dataSeries.Count() == 0 || chart == null || chart.Series == null)
            {
                return;
            }
            if (chart.Series != null && chart.Series.Count == 0)
            {
                chart.Series.Add(new PieSeries());
            }
            ((PieSeries)chart.Series[0]).ItemsSource = dataSeries;
        }


        /// <summary>
        /// Format piechart data (converts to two values - one for in-ploidy and one for out-of-ploidu
        /// </summary>
        private static KeyValuePair<string, double>[] GetPieData(Collection<double> data, int ploidyLevel)
        {
            if (data == null || data.Count() == 0)
            {
                return null;
            }
            KeyValuePair<string, double>[] dataSeries = new KeyValuePair<string, double>[2];
            int j = 0;
            double inPloidy = 0;
            double outOfPloidy = 0;
            foreach (double d in data)
            {
                if (j < ploidyLevel)
                {
                    inPloidy += d;
                }
                else if (j >= ploidyLevel)
                {
                    outOfPloidy += d;
                }

                if (j == ploidyLevel - 1)
                {
                    dataSeries[0] = new KeyValuePair<string, double>("In ploidy", inPloidy);
                }
                else if (j == data.Count() - 1)
                {
                    dataSeries[1] = new KeyValuePair<string, double>("Out of ploidy", outOfPloidy);
                }
                j++;
            }
            return dataSeries;
        }

     
        /// <summary>
        /// Update SliderB and columnchart x axis when SliderA value changes
        /// </summary>
        private void SliderA_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            ReadCountX.Minimum = e.NewValue;
            SliderB.Minimum = e.NewValue + 1;
        }

        /// <summary>
        /// Update SliderA and columnchart x axis when SliderB value changes
        /// </summary>
        private void SliderB_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            ReadCountX.Maximum = e.NewValue;
            SliderA.Maximum = e.NewValue - 1;
        }

        /// <summary>
        /// Get a list of all charts which are displayed on-screen
        /// </summary>
        /// <returns></returns>
        private List<Chart> GetDisplayedChart()
        {
            List<Chart> charts = new List<Chart>();
            if (Tab1.IsSelected)
            {
                charts.Add(PloidyPie);
                charts.Add(PloidyPieGood);
            }
            else if (Tab2.IsSelected)
            {
                charts.Add(ReadCountChart);
            }
            else if (Tab3.IsSelected)
            {
                charts.Add(IndivChart);
            }
            else
            {
                return null;
            }
            return charts;
        }

        #endregion

        #region process sequences

        /// <summary>
        /// Initialise handler and handler settings based on user's input
        /// </summary>
        private void InitHandler(string filename, string ploidy, string dirtCutoff, string ploidyCutoff, string alignQualCutoff,
            string readQualCutoff, string popPercent, string hapMaxCutoff, bool? onlyHaplotypeGood, string numSamples, bool? outputToFile, bool? metricFileParent,
            bool? metricFileChild, bool? outputOverviewParent, bool? outputOverviewChild, bool? genotypesToFile, bool? haplotypesToFile)
        {
            // Check user input
            string newName = filename.Split(new char[] { '.' })[0];

            // Set up the handler
            metricHandler = new ClusterMetricHandler(newName, Convert.ToInt32(ploidy, ci), Convert.ToInt32(numSamples, ci));

            metricHandler.DirtCutoff = Convert.ToDouble(dirtCutoff, ci);
            metricHandler.AlignmentQualityCutoff = Convert.ToDouble(alignQualCutoff, ci);
            metricHandler.ReadQualityCutoff = Convert.ToDouble(readQualCutoff, ci);
            metricHandler.PopulationPercentageCutoff = Convert.ToDouble(popPercent, ci);

            metricHandler.PloidyDisagreementCutoff = Convert.ToDouble(ploidyCutoff, ci);
            metricHandler.HaplotypesMaxCutoff = Convert.ToInt32(hapMaxCutoff, ci);
            metricHandler.OnlyHaplotypeGood = onlyHaplotypeGood == true;
            metricHandler.WriteToFilteredBam = outputToFile == true;
            metricHandler.WriteClusterMetricOriginal = metricFileParent == true;
            metricHandler.WriteClusterMetricFiltered = metricFileChild == true;
            metricHandler.WriteOverviewMetricOriginal = outputOverviewParent == true;
            metricHandler.WriteOverviewMetricFiltered = outputOverviewChild == true;
            metricHandler.WriteGenotypesFile = genotypesToFile == true;
            metricHandler.WriteHaplotypesFile = haplotypesToFile == true;
        }

        /// <summary>
        /// Read input file header and set progress bar based on number of sequences in input file
        /// </summary>
        private void ReadHeader(string filename, BAMParser parser, int numClustersInInputFile)
        {
            using (Stream readStream = new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                metricHandler.InputHeader = parser.GetHeader(readStream);
                numClustersInInputFile = metricHandler.InputHeader.ReferenceSequences.Count;
                Console.WriteLine(Properties.Resources.CLUSTER_COUNT_DISPLAY + numClustersInInputFile);
            }
            Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
                new IntDelegate(SetProgressBar), numClustersInInputFile);
        }

        /// <summary>
        /// Create handler and initialise handler settings, then process sequences
        /// </summary>
        private System.Delegate ParseBAMMetric(string filename, string ploidy, string dirtCutoff, string ploidyDisagreement,
            string alignQualCutoff, 
            string readQualCutoff, string popPercent, string hapMaxCutoff, bool? onlyHaplotypeGood, string numSamples, bool? outputToFile, 
            bool? metricFileParent, bool? metricFileChild, bool? outputOverviewParent, bool? outputOverviewChild, 
            bool? genotypesToFile, bool? haplotypesToFile)
        {
            using (BAMParser parser = new BAMParser())
            {
                int numClustersInInputFile = 0;
            
                // Initialise the metric handler
                InitHandler(filename, ploidy, dirtCutoff, ploidyDisagreement, alignQualCutoff, readQualCutoff, popPercent, hapMaxCutoff, onlyHaplotypeGood, numSamples, 
                    outputToFile, metricFileParent, metricFileChild, outputOverviewParent, outputOverviewChild, genotypesToFile, haplotypesToFile);

                ReadHeader(filename, parser, numClustersInInputFile);

                // Begin parsing sequences
                Dispatcher.BeginInvoke( System.Windows.Threading.DispatcherPriority.Normal,
                    new QuickDelegate(UpdateGui_BeganParsing)); // updates GUI to indicate that parsing has begun
                ProcessSequences(filename, parser);
                
                // Finished sequences for this input file
                Dispatcher.BeginInvoke( System.Windows.Threading.DispatcherPriority.Normal,
                    new QuickDelegate(UpdateGui_FinishedParsing));

                // Close handler and return
                Console.WriteLine(Properties.Resources.FINISHED);
                metricHandler.Dispose();
            }
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
                if (metricHandler.Add(se))
                {
                    clusterCount = metricHandler.ClusterCount;
                    increment = (clusterCount < increaseIncrementAfter) ? incrementOne : incrementTwo;
                    if ((clusterCount == 1 || (double)clusterCount % increment == 0) && clusterCount != updateDisplayForClusterIndex)
                    {
                        updateDisplayForClusterIndex = clusterCount;

                        Dispatcher.BeginInvoke(
                            System.Windows.Threading.DispatcherPriority.Normal,
                            new StatsDelegate(UpdateStatsPanel), metricHandler.GoodCount, clusterCount,
                                metricHandler.MaxSampleCount, metricHandler.MaxAlignmentQuality, metricHandler.MaxReadQuality,
                                metricHandler.AverageDirt, metricHandler.AverageMapQ, metricHandler.AverageReadQ,
                                metricHandler.AverageDirtGood, metricHandler.AverageMapQGood, metricHandler.AverageReadQGood);
                    }
                }
                else
                {
                    Console.WriteLine(Properties.Resources.HANDLER_FINISHED);
                    break;
                }
            }
            // Tell the handler that there are no more sequences to receive
            metricHandler.SetComplete(); 
        }

        #endregion

        #region menu

        /// <summary>
        /// Display BEGIN menu
        /// </summary>
        private void ShowBeginMenu(object sender, RoutedEventArgs e)
        {
            BeginMenu.Visibility = System.Windows.Visibility.Visible;
            SaveMenu.Visibility = System.Windows.Visibility.Hidden;
            HelpMenu.Visibility = System.Windows.Visibility.Hidden;
            AbortMenu.Visibility = System.Windows.Visibility.Hidden;
        }

        /// <summary>
        /// Display HELP menu
        /// </summary>
        private void ShowHelpMenu(object sender, RoutedEventArgs e)
        {
            BeginMenu.Visibility = System.Windows.Visibility.Hidden;
            SaveMenu.Visibility = System.Windows.Visibility.Hidden;
            HelpMenu.Visibility = System.Windows.Visibility.Visible;
            AbortMenu.Visibility = System.Windows.Visibility.Hidden;
        }

        /// <summary>
        /// Display SAVE menu
        /// </summary>
        private void ShowSaveMenu(object sender, RoutedEventArgs e)
        {
            BeginMenu.Visibility = System.Windows.Visibility.Hidden;
            SaveMenu.Visibility = System.Windows.Visibility.Visible;
            HelpMenu.Visibility = System.Windows.Visibility.Hidden;
            AbortMenu.Visibility = System.Windows.Visibility.Hidden;
        }

        /// <summary>
        /// Display ABORT menu
        /// </summary>
        private void ShowAbortMenu(object sender, RoutedEventArgs e)
        {
            BeginMenu.Visibility = System.Windows.Visibility.Hidden;
            SaveMenu.Visibility = System.Windows.Visibility.Hidden;
            HelpMenu.Visibility = System.Windows.Visibility.Hidden;
            AbortMenu.Visibility = System.Windows.Visibility.Visible;
        }

        /// <summary>
        /// Close BEGIN menu
        /// </summary>
        private void Begin_Menu_Close(object sender, RoutedEventArgs e)
        {
            BeginMenu.Visibility = System.Windows.Visibility.Hidden;
        }
        
        /// <summary>
        /// Close ABORT menu
        /// </summary>
        private void Abort_Menu_Close(object sender, RoutedEventArgs e)
        {
            AbortMenu.Visibility = System.Windows.Visibility.Hidden;
        }
        
        /// <summary>
        /// Close HELP menu
        /// </summary>
        private void Help_Menu_Close(object sender, RoutedEventArgs e)
        {
            HelpMenu.Visibility = System.Windows.Visibility.Hidden;
        }
        
        /// <summary>
        /// Close SAVE menu
        /// </summary>
        private void Save_Menu_Close(object sender, RoutedEventArgs e)
        {
            SaveMenu.Visibility = System.Windows.Visibility.Hidden;
        }
        
        /// <summary>
        /// Show INPUT menu
        /// </summary>
        private void ShowInputMenu(object sender, RoutedEventArgs e)
        {
            Tab0.IsSelected = true;
            FileExpander.IsExpanded = true;
            Begin_Menu_Close(null,null);
            FileExpander.Focus();
        }

        /// <summary>
        /// Show OUTPUT menu
        /// </summary>
        private void ShowOutputMenu(object sender, RoutedEventArgs e)
        {
            Tab0.IsSelected = true;
            OutputFileExpander.IsExpanded = true;
            Begin_Menu_Close(null,null);
            OutputFileExpander.Focus();
        }

        /// <summary>
        /// Show HELP topic: Input
        /// </summary>
        private void ShowHelpInput(object sender, RoutedEventArgs e)
        {
            Tab4.IsSelected = true;
            AboutInputParameters.IsExpanded = true;
            Help_Menu_Close(null,null);
            AboutInputParameters.Focus();
        }

        /// <summary>
        /// Show HELP topic: Output
        /// </summary>
        private void ShowHelpOutput(object sender, RoutedEventArgs e)
        {
            Tab4.IsSelected = true;
            AboutOutputFiles.IsExpanded = true;
            Help_Menu_Close(null,null);
            AboutOutputFiles.Focus();
        }

        /// <summary>
        /// Show HELP topic: About
        /// </summary>
        private void ShowHelpAbout(object sender, RoutedEventArgs e)
        {
            Tab4.IsSelected = true;
            AboutAbout.IsExpanded = true;
            Help_Menu_Close(null,null);
            AboutAbout.Focus();
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
                MessageBoxResult result = MessageBox.Show(Properties.Resources.WINDOW_CLOSE_WARNING,
                Properties.Resources.GENERAL_WARNING, MessageBoxButton.YesNo, MessageBoxImage.Warning);
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
            TextBox t = sender as TextBox;
            string s = "";
            int dec = 0;
            foreach (Char ch in t.Text.ToCharArray())
            {
                s += Char.IsDigit(ch) || (ch == '.' && dec == 0)? ch.ToString() : "";
                if (ch == '.')
                {
                    dec++;
                }
            }
            t.Text = (String.IsNullOrEmpty(s)) ? "0" : s; // if blank, set to 0
        }

        #endregion

        #region bind labels to input fields

        /// <summary>
        /// Bind click on label to click on checkbox
        /// </summary>
        private void Label_MouseDown(object sender, MouseButtonEventArgs e)
        {
            OutputToFile.IsChecked = (OutputToFile.IsEnabled) ? !OutputToFile.IsChecked : OutputToFile.IsChecked;
        }

        /// <summary>
        /// Bind click on label to click on checkbox
        /// </summary>
        private void Label_MouseDown_1(object sender, MouseButtonEventArgs e)
        {
            OutputOverviewChild.IsChecked = (OutputOverviewChild.IsEnabled) ? !OutputOverviewChild.IsChecked : OutputOverviewChild.IsChecked;
        }

        /// <summary>
        /// Bind click on label to click on checkbox
        /// </summary>
        private void Label_MouseDown_2(object sender, MouseButtonEventArgs e)
        {
            OutputOverviewParent.IsChecked = (OutputOverviewParent.IsEnabled) ? !OutputOverviewParent.IsChecked : OutputOverviewParent.IsChecked;
        }

        /// <summary>
        /// Bind click on label to click on checkbox
        /// </summary>
        private void Label_MouseDown_3(object sender, MouseButtonEventArgs e)
        {
            MetricFileChild.IsChecked = (MetricFileChild.IsEnabled) ? !MetricFileChild.IsChecked : MetricFileChild.IsChecked;
        }

        /// <summary>
        /// Bind click on label to click on checkbox
        /// </summary>
        private void Label_MouseDown_4(object sender, MouseButtonEventArgs e)
        {
            MetricFileParent.IsChecked = (MetricFileParent.IsEnabled) ? !MetricFileParent.IsChecked : MetricFileParent.IsChecked;
        }

        /// <summary>
        /// Bind click on label to click on checkbox
        /// </summary>
        private void Label_MouseDown_5(object sender, MouseButtonEventArgs e)
        {
            GenotypesToFile.IsChecked = (GenotypesToFile.IsEnabled) ? !GenotypesToFile.IsChecked : GenotypesToFile.IsChecked;
        }

        /// <summary>
        /// Bind click on label to click on checkbox
        /// </summary>
        private void Label_MouseDown_6(object sender, MouseButtonEventArgs e)
        {
            HaplotypesToFile.IsChecked = (HaplotypesToFile.IsEnabled) ? !HaplotypesToFile.IsChecked : HaplotypesToFile.IsChecked;
        }

        /// <summary>
        /// Bind click on label to click on checkbox
        /// </summary>
        private void Label_MouseDown_7(object sender, MouseButtonEventArgs e)
        {
            OnlyHaplotypeGood.IsChecked = (OnlyHaplotypeGood.IsEnabled) ? !OnlyHaplotypeGood.IsChecked : OnlyHaplotypeGood.IsChecked;
        }

        #endregion

        #region extra features

        /// <summary>
        /// Save displayed chart as an image
        /// </summary>
        private void SaveChartImage(object sender, RoutedEventArgs e)
        {
            List<Chart> charts = GetDisplayedChart();
            if (charts != null)
            {
                foreach (Chart chart in charts)
                {
                    RenderTargetBitmap renderBitmap = new RenderTargetBitmap(
                         (int)chart.ActualWidth,
                         (int)chart.ActualHeight,
                         96d,
                         96d,
                         PixelFormats.Pbgra32);

                    renderBitmap.Render(chart);
                    string path = metricHandler.FileName + "\\" + chart.Name + ".bmp";
                    using (FileStream outStream = new FileStream(path, FileMode.Create))
                    {
                        PngBitmapEncoder encoder = new PngBitmapEncoder();
                        encoder.Frames.Add(BitmapFrame.Create(renderBitmap));
                        encoder.Save(outStream);
                    }
                }
            }
        }

        #endregion

        private void ExpectedPloidyTextbox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (HaplotypesToFile != null) // only perform if the rest of the form has been initialised
            {
                ComboBox s = (ComboBox)sender;
                if(s.SelectedIndex == 1) // diploid
                {
                    MaxHaplotypesTextbox.IsEnabled = true;
                    OnlyHaplotypeGood.IsEnabled = true;
                    GenotypesToFile.IsEnabled = true;
                    HaplotypesToFile.IsEnabled = true;
                }
                else
                {
                    MaxHaplotypesTextbox.IsEnabled = false;
                    OnlyHaplotypeGood.IsEnabled = false;
                    GenotypesToFile.IsEnabled = false;
                    HaplotypesToFile.IsEnabled = false;

                    GenotypesToFile.IsChecked = false;
                    HaplotypesToFile.IsChecked = false;
                }
            }
        }

        #endregion

    }
}
