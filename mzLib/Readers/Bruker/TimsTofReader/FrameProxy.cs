﻿using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using ThermoFisher.CommonCore.Data.Interfaces;

namespace Readers.Bruker.TimsTofReader
{
    internal class FrameProxyFactory
    {
        internal FrameTable FramesTable { get; }
        internal UInt64 FileHandle { get; }
        internal Object FileLock { get; }
        internal TimsConversion Converter { get; }

        internal FrameProxyFactory(FrameTable table, UInt64 fileHandle, Object fileLock)
        {
            FramesTable = table;
            FileHandle = fileHandle;
            FileLock = fileLock;
            Converter = new TimsConversion(fileHandle, fileLock);
        }   

        internal FrameProxy GetFrameProxy(long frameId)
        {
            return new FrameProxy(FileHandle, frameId, FramesTable.NumScans[frameId - 1], FileLock, Converter);
        }

        internal Polarity GetPolarity(long frameId)
        {
            return FramesTable.Polarity[frameId - 1] == '+' ? Polarity.Positive : Polarity.Negative;
        }

        internal double GetRetentionTime(long frameId)
        {
            return (double)FramesTable.RetentionTime[frameId - 1];
        }

        internal double GetInjectionTime(long frameId)
        {
            return FramesTable.FillTime[frameId - 1];
        }

        internal int GetScanCount(long frameId)
        {
            return FramesTable.NumScans[frameId - 1];
        }

        internal double GetInjectionTimeSum(long firstFrameId, long lastFrameId)
        {
            double injectionTimeSum = 0;
            for(long i = firstFrameId; i <= lastFrameId; i++)
            {
                injectionTimeSum += FramesTable.FillTime[i - 1];
            }
            return injectionTimeSum;
        }
    }

    [StructLayout(LayoutKind.Sequential)]
    public unsafe struct CentroidedSpectrum
    {

        internal IntPtr MzPointer { get; set; }
        internal IntPtr AreaPointer { get; set; }
        internal UInt32 NumPeaks { get; set; }

        public CentroidedSpectrum()
        {
            MzPointer = IntPtr.Zero;
            AreaPointer = IntPtr.Zero;
            NumPeaks = 0;
        }

    }

    public static unsafe class MsmsSpectrumWrapper
    {
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        public delegate void msms_spectrum_function(Int64 id, UInt32 num_peaks, in double*[] mz_values, in float*[] area_values, void* user_data);

        /// <summary>
        /// Function type that takes a centroided peak list.
        /// </summary>
        /// <param name="id"> the id of the precursor or frame </param>
        /// <param name="num_peaks"> the number of peaks in the spectrum </param>
        /// <param name="mz_values"> all peak m/z values </param>
        /// <param name="area_values"> all peak areas </param>
        /// <param name="user_data"> Passed to the function by "tims_read_pasef_msms_v2" function </param>
        internal unsafe static void GetCentroidedSpectrumCallback(Int64 id, UInt32 num_peaks, in double*[] mz_values, in float*[] area_values, void* user_data)
        {
            //int mz_buffer_size = (int)num_peaks * Marshal.SizeOf<double>();
            //int area_buffer_size = (int)num_peaks * Marshal.SizeOf<float>();
            //IntPtr mz_pointer = Marshal.AllocHGlobal(mz_buffer_size);
            //IntPtr area_pointer = Marshal.AllocHGlobal(area_buffer_size);

            //var mzArray = new double[mz_buffer_size];
            //FrameProxy.CopyToManaged(mz_pointer, mzArray, 0, mz_buffer_size);
            //Marshal.FreeHGlobal(mz_pointer);

            //var areaArray = new float[area_buffer_size];
            //FrameProxy.CopyToManaged(area_pointer, areaArray, 0, area_buffer_size);
            //Marshal.FreeHGlobal(area_pointer);

            //CentroidedSpectrum spectrum = Marshal.PtrToStructure<CentroidedSpectrum>((IntPtr)user_data);
            //spectrum.MzPointer = mz_pointer;
            //spectrum.AreaPointer = area_pointer;
            //spectrum.NumPeaks = num_peaks;
        }

        [MarshalAs(UnmanagedType.FunctionPtr)]
        private static msms_spectrum_function? callback;

        public static unsafe UInt32 GetSpectrumTest(UInt64 handle, Int64 frame_id, UInt32 scan_begin, UInt32 scan_end)
        {
            CentroidedSpectrum spectrum = new CentroidedSpectrum();
            IntPtr spectrumPointer = Marshal.AllocHGlobal(Marshal.SizeOf<CentroidedSpectrum>());
            Marshal.StructureToPtr(spectrum, spectrumPointer, false);

            callback = new msms_spectrum_function(GetCentroidedSpectrumCallback);
            //IntPtr callbackPtr = GCHandle.ToIntPtr(GCHandle.Alloc(callback));
            var fPtr = Marshal.GetFunctionPointerForDelegate(callback);


            UInt32 returnCode = tims_extract_centroided_spectrum_for_frame_v2(handle, frame_id, scan_begin, scan_end, fPtr, spectrumPointer);

            CentroidedSpectrum testSpec = Marshal.PtrToStructure<CentroidedSpectrum>(spectrumPointer);

            return returnCode;
        }


        public static MzSpectrum GetMzSpectrum(CentroidedSpectrum spectrum)
        {
            int mz_buffer_size = (int)spectrum.NumPeaks * Marshal.SizeOf<double>();
            int area_buffer_size = (int)spectrum.NumPeaks * Marshal.SizeOf<float>();
            IntPtr mz_pointer = Marshal.AllocHGlobal(mz_buffer_size);
            IntPtr area_pointer = Marshal.AllocHGlobal(area_buffer_size);

            var mzArray = new double[mz_buffer_size];
            FrameProxy.CopyToManaged(mz_pointer, mzArray, 0, mz_buffer_size);
            Marshal.FreeHGlobal(mz_pointer);

            var areaArray = new float[area_buffer_size];
            FrameProxy.CopyToManaged(area_pointer, areaArray, 0, area_buffer_size);
            Marshal.FreeHGlobal(area_pointer);

            var areaDoubleArray = Array.ConvertAll(areaArray, entry => (double)entry);

            return new MzSpectrum(mzArray, areaDoubleArray, shouldCopy: false);
        }

        [DllImport("Bruker/TimsTofReader/timsdata.dll", CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern UInt32 tims_extract_centroided_spectrum_for_frame_v2
              (UInt64 handle, Int64 frame_id, UInt32 scan_begin, UInt32 scan_end, IntPtr callback, IntPtr user_data);
    }

    internal class FrameProxy
    {
        private int[] _scanOffsets; // Number of peaks that precede a given scan in a frame
        private uint[] _rawData;
        private const int _defaultBufferSize = 4096;
        internal UInt64 FileHandle { get; }
        internal long FrameId { get; }
        internal int NumberOfScans { get; }
        internal int TotalNumberOfPeaks => _scanOffsets[NumberOfScans - 1];
        internal TimsConversion Converter { get; }

        internal FrameProxy(UInt64 fileHandle, long frameId, int numScans, Object fileLock, TimsConversion converter)
        {
            NumberOfScans = numScans;
            FileHandle = fileHandle;
            FrameId = frameId;
            Converter = converter;
            _rawData = GetScanRawData(fileHandle, frameId, (uint)numScans, fileLock);
            _scanOffsets = PartialSum(_rawData, 0, numScans);
        }

        internal double[] GetScanMzs(int zeroIndexedScanNumber)
        {
            double[] mzLookupIndices = Array
                .ConvertAll(_rawData[GetXRange(zeroIndexedScanNumber)], entry => (double)entry);

            return Converter.DoTransformation(FileHandle, FrameId, mzLookupIndices, ConversionFunctions.IndexToMz);
        }

        internal int[] GetScanIntensities(int zeroIndexedScanNumber)
        {
            return Array.ConvertAll(_rawData[GetYRange(zeroIndexedScanNumber)], entry => (int)entry);
        }

        internal double GetOneOverK0(double zeroIndexedScanNumberMedian)
        {
            ThrowIfInvalidScanNumber((int)zeroIndexedScanNumberMedian);
            double[] scanNumberArray = new double[] { zeroIndexedScanNumberMedian };
            return Converter
                .DoTransformation(FileHandle, FrameId, scanNumberArray, ConversionFunctions.ScanToOneOverK0)[0];
        }

        /// <summary>
        /// Read a range of scans from a single frame.
        ///
        /// Output layout: (N = scan_end - scan_begin = number of requested scans)
        ///   N x uint32_t: number of peaks in each of the N requested scans
        ///   N x (two uint32_t arrays: first indices, then intensities)
        ///
        /// Note: different threads must not read scans from the same storage handle
        /// concurrently.
        /// </summary> 
        internal unsafe static uint[] GetScanRawData(UInt64 fileHandle, long frameId, UInt32 numScans, Object fileLock)
        {
            int bufferSize = _defaultBufferSize;
            // buffer expansion loop
            while (true)
            {
                IntPtr pData = Marshal.AllocHGlobal(bufferSize * Marshal.SizeOf<Int32>());
                uint outputLength;

                lock ( fileLock )
                {
                    outputLength = tims_read_scans_v2(
                    fileHandle,
                    frameId,
                    scan_begin: 0,
                    scan_end: numScans,
                    buffer: pData,
                    length: (uint)(bufferSize * 4));
                }
                
                if (4 * bufferSize > outputLength)
                {
                    var dataArray = new uint[bufferSize];
                    CopyToManaged(pData, dataArray, 0, bufferSize);
                    Marshal.FreeHGlobal(pData);

                    return dataArray;
                }

                if (outputLength > 16777216) // Arbitrary 16 mb frame limit
                {
                    throw new Exception("Maximum frame size exceeded");
                }

                // Increase buffer size if necessary
                bufferSize = ((int)outputLength / 4) + 1;

                Marshal.FreeHGlobal(pData);
            }
        }

        internal int GetNumberOfPeaks(int zeroIndexedScanNumber)
        {
            ThrowIfInvalidScanNumber(zeroIndexedScanNumber);
            return (int)_rawData[zeroIndexedScanNumber];
        }

        /// <summary>
        /// Returns a range containing the start(inclusive) and end (exclusive) indices
        /// for the segment of the _rawData array corresponding to the m/z lookup values for 
        /// a given scan
        /// </summary>
        /// <exception cref="ArgumentException"> Throws exception if scan number out of range </exception>
        internal Range GetXRange(int zeroIndexedScanNumber)
        {
            ThrowIfInvalidScanNumber(zeroIndexedScanNumber);
            return GetScanRange(zeroIndexedScanNumber, offset: 0);
        }

        /// <summary>
        /// Returns a range containing the start(inclusive) and end (exclusive) indices
        /// for the segment of the _rawData array corresponding to raw intensity values for a given scan
        /// </summary>
        internal Range GetYRange(int zeroIndexedScanNumber)
        {
            ThrowIfInvalidScanNumber(zeroIndexedScanNumber);
            return GetScanRange(zeroIndexedScanNumber, offset: (int)_rawData[zeroIndexedScanNumber]);
        }

        /// <exception cref="ArgumentException"> Throws exception if scan number out of range </exception>
        private void ThrowIfInvalidScanNumber(int zeroIndexedScanNumber)
        {
            if (zeroIndexedScanNumber < 0 || zeroIndexedScanNumber >= NumberOfScans)
                throw new ArgumentException("Scan number out of range.");
        }

        private Range GetScanRange(int zeroIndexedScanNumber, int offset)
        {
            int start = NumberOfScans + 2*_scanOffsets[zeroIndexedScanNumber] + offset;
            return new Range(start, start + (int)_rawData[zeroIndexedScanNumber]);
        }

        /// <summary>
        /// Calculates the running total of an array, beginning with 
        /// the start index (inclusive) and ending with the end index (exclusive).
        /// Used for determining scan offsets.
        /// </summary>
        /// <param name="array">Array to be summed</param>
        /// <param name="start">Where to begin summing</param>
        /// <param name="end">Where summing ends (exclusive) </param>
        /// <returns> An array of length (end - start) containing the 
        /// partial sums at each index of the input array</returns>
        public static int[] PartialSum(uint[] array, int start, int end)
        {
            int runningTotal = 0;
            int[] sums = new int[end - start + 1];
            sums[0] = 0;

            for(int i = 0; i < end; i++)
            {
                runningTotal += (int)array[i];
                sums[i+1] = runningTotal;
            }
            return sums;
        }

        /// <summary>
        /// This is reimplementation of the Marshal.Copy method that allows for arbitrary types
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="source"></param>
        /// <param name="destination"></param>
        /// <param name="startIndex"></param>
        /// <param name="length"></param>
        /// <exception cref="ArgumentNullException"></exception>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        internal static unsafe void CopyToManaged<T>(IntPtr source, T[] destination, int startIndex, int length)
        {
            if (source == IntPtr.Zero) throw new ArgumentNullException(nameof(source));
            if (destination is null) throw new ArgumentNullException(nameof(destination));
            if (startIndex < 0) throw new ArgumentOutOfRangeException(nameof(startIndex));
            if (length < 0) throw new ArgumentOutOfRangeException(nameof(length));

            void* sourcePtr = (void*)source;
            Span<T> srcSpan = new Span<T>(sourcePtr, length);
            Span<T> destSpan = new Span<T>(destination, startIndex, length);

            srcSpan.CopyTo(destSpan);
        }


        /// <summary>
        /// Read a range of scans from a single frame.
        ///
        /// Output layout: (N = scan_end - scan_begin = number of requested scans)
        ///   N x uint32_t: number of peaks in each of the N requested scans
        ///   N x (two uint32_t arrays: first indices, then intensities)
        ///
        /// Note: different threads must not read scans from the same storage handle
        /// concurrently.
        /// </summary> 
        /// <param name="handle"> Unique Handle of .d file ( returned on tims_open() )</param>
        /// <param name="frame_id"> From .tdf SQLite: Frames.Id </param>
        /// <param name="scan_begin"> first scan number to read (inclusive) </param>
        /// <param name="scan_end"> Last scan number (exclusive) </param>
        /// <param name="buffer"> Destination buffer allocated by user </param>
        /// <param name="length"> Length of the buffer (in bytes, i.e. 4 * buffer.length) </param>
        /// <returns> 0 on error, otherwise the number of buffer bytes necessary for the output
        /// of this call (if this is larger than the provided buffer length, the result is not
        /// complete). </returns>
        [DllImport("Bruker/TimsTofReader/timsdata.dll", CallingConvention = CallingConvention.Cdecl)]
        unsafe static extern UInt32 tims_read_scans_v2
              (UInt64 handle, Int64 frame_id, UInt32 scan_begin, UInt32 scan_end, IntPtr buffer, UInt32 length);

    }
}