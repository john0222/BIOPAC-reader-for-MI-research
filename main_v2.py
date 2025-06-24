import bioread
import numpy as np
from datetime import timedelta
from tpalgo import RatECGAnalyzer
import matplotlib.pyplot as plt

from pylab import mpl
mpl.rcParams["font.sans-serif"] = ["SimHei"]
mpl.rcParams["axes.unicode_minus"] = False

def read_acq_data(filepath):
    """
    read from acq file 

    """
    try:
        data = bioread.read_file(filepath)
        print(f"成功讀取檔案: {filepath}")
        print(f"包含 {len(data.channels)} 個通道。")
        return data
    except Exception as e:
        print(f"讀取檔案時發生錯誤: {e}")
        return None

def get_channel_data(data_manager, channel_name_keywords):
    """
    Finds a channel by keywords, prioritizing exact matches.

    Args:
        data_manager (bioread.read_file): The data manager object.
        channel_name_keywords (list): A list of keywords to search for.

    Returns:
        bioread.Channel: The channel object if found, None otherwise.
    """
    for ch in data_manager.channels:
        if ch.name.lower() in [kw.lower() for kw in channel_name_keywords]:
            print(f"找到通道 (精確匹配): {ch.name}")
            return ch
    for ch in data_manager.channels:
        for kw in channel_name_keywords:
            if kw.lower() in ch.name.lower():
                print(f"找到通道 (部分匹配): {ch.name}")
                return ch
    print(f"未找到符合關鍵字 {channel_name_keywords} 的通道。")
    return None

def calculate_blood_pressure_metrics(bp_channel, sample_rate, start_second, time_intervals):
    """
    Calculates blood pressure metrics (instantaneous max/min, average) at specified time points.

    Args:
        bp_channel (bioread.Channel): The blood pressure channel object.
        sample_rate (float): The sampling rate of the channel.
        start_second (int): The starting point in seconds for calculations.
        time_intervals (list): A list of time intervals (in seconds) from the start_second
                                at which to calculate metrics.

    Returns:
        list: A list of dictionaries, each containing 'time_point', 'instantaneous_bp_max',
              'instantaneous_bp_min', and 'average_bp' for each interval.
    """
    bp_data = bp_channel.data
    results = []

    print(f"\n計算血壓指標 (從 {start_second} 秒開始)...")

    # Define a small window for instantaneous BP (e.g., 0.5 seconds)
    instantaneous_window_samples = int(0.5 * sample_rate)

    for interval_seconds in time_intervals:
        current_time_seconds = start_second + interval_seconds
        current_sample_index = int(current_time_seconds * sample_rate)

        # Ensure index is within data bounds
        if current_sample_index >= len(bp_data):
            print(f"警告: 時間點 {current_time_seconds} 秒超出數據範圍，跳過。")
            continue

        # Calculate instantaneous BP (max/min) within a small window around the current_sample_index
        start_instantaneous_idx = max(0, current_sample_index - instantaneous_window_samples // 2)
        end_instantaneous_idx = min(len(bp_data), current_sample_index + instantaneous_window_samples // 2)

        instantaneous_bp_segment = bp_data[start_instantaneous_idx:end_instantaneous_idx]
        instantaneous_bp_max = np.max(instantaneous_bp_segment) if instantaneous_bp_segment.size > 0 else np.nan
        instantaneous_bp_min = np.min(instantaneous_bp_segment) if instantaneous_bp_segment.size > 0 else np.nan

        # Calculate average BP from the start_second up to the current_time_seconds
        average_bp_segment = bp_data[int(start_second * sample_rate):current_sample_index]
        average_bp = np.mean(average_bp_segment) if average_bp_segment.size > 0 else np.nan

        results.append({
            'time_point': timedelta(seconds=current_time_seconds),
            'instantaneous_bp_max': instantaneous_bp_max,
            'instantaneous_bp_min': instantaneous_bp_min,
            'average_bp': average_bp
        })
    return results
#--------------------------------------------------------------------------------------------------------------------------
def analyze_ecg_for_arrhythmias(file ,ecg_channel, sample_rate, start_second, time_intervals):
  
    
    ecg_data = ecg_channel.data
    all_ecg_results = []

#---------rr analyzer
    acq_file_path = file
    try:

        data = bioread.read_file(acq_file_path)


        if len(data.channels) > 0:
            #search for the channel
            ecg_channel = None
            for channel in data.channels:
                print("search ing for data channel")
                if 'ecg' in channel.name.lower() or 'heart' in channel.name.lower():
                    ecg_channel = channel
                    break
            
            if ecg_channel is None and len(data.channels) > 0:
                print("未找到明確的 'ECG' 或 'Heart' 通道，將使用第一個通道。")
                ecg_channel = data.channels[0] # default on first channel 

            if ecg_channel:
                ecg_signal = ecg_channel.data.flatten()
                sfreq = ecg_channel.samples_per_second # sample rate 

                print(f"成功讀取通道: {ecg_channel.name}")
                print(f"採樣頻率: {sfreq} Hz")
                print(f"ECG 訊號長度: {len(ecg_signal)} 點")

                # create RatECGAnalyzer object
                # 根據大鼠的ECG特性調整參數
                analyzer = RatECGAnalyzer(
                    sfreq=sfreq,
                    highpass_cutoff=8.0,       # 提高高通截止頻率
                    lowpass_cutoff=50.0,       # 提高低通截止頻率
                    integration_window_ms=80,  # 縮短集成窗口
                    threshold_factor=0.3,
                    noise_threshold_factor=0.6,
                    refractory_period_ms=100
                )

                # EKG analyze
                r_peak_indices, rr_intervals_ms = analyzer.analyze_ecg(ecg_signal)

                print(f"\n檢測到的R波數量: {len(r_peak_indices)}")
                if len(r_peak_indices) > 0:
                    print(f"檢測到的R波索引: {r_peak_indices[:10]}...")
                if len(rr_intervals_ms) > 0:
                    print(f"計算的R-R間期 (毫秒): {rr_intervals_ms[:10]}...")

                # plotting
                plt.figure(figsize=(15, 6))
                time_axis = np.arange(len(ecg_signal)) / sfreq
                plt.plot(time_axis, ecg_signal, label='原始ECG訊號', alpha=0.7)
                plt.plot(time_axis[r_peak_indices], ecg_signal[r_peak_indices], 'ro', markersize=6, label='檢測到的R波')
                plt.title(f'ECG R-R test result on channel: {ecg_channel.name})')
                plt.xlabel('time(sec)')
                plt.ylabel('amplititude')
                plt.grid(True)
                plt.legend()
                plt.tight_layout()
                plt.show()

                # RR interval plot
                if rr_intervals_ms:
                    plt.figure(figsize=(10, 5))
                    plt.hist(rr_intervals_ms, bins=20, edgecolor='black')
                    plt.title('R-R interval spread')
                    plt.xlabel('R-R interval length')
                    plt.ylabel('frequency')
                    plt.grid(True)
                    plt.tight_layout()
                    plt.show()

                    # calculate average bpm
                    average_rr_ms = np.mean(rr_intervals_ms)
                    average_hr_bpm = 60000 / average_rr_ms if average_rr_ms > 0 else 0
                    print(f"average R-R interval: {average_rr_ms:.2f} ms")
                    print(f"average BPM: {average_hr_bpm:.2f} BPM")
                else:
                    print("NOT ENOUGH R-R INTERVAL DATA FOR CALCULATION")

            else:
                print("no ECG channel found in file")
        else:
            print(f"ERROR : file '{acq_file_path}' no channal found")

    except FileNotFoundError:
        print(f"錯誤: 找不到檔案 '{acq_file_path}'。請檢查路徑是否正確。")
    except Exception as e:
        print(f"讀取或處理檔案時發生錯誤: {e}")

#---------------------


    print(f"\n分析ECG心律不整 (從 {start_second} 秒開始，待完成)...")

    for interval_seconds in time_intervals:
        current_time_seconds = start_second + interval_seconds
        current_sample_index = int(current_time_seconds * sample_rate)

        if current_sample_index >= len(ecg_data):
            print(f"警告: ECG時間點 {current_time_seconds} 秒超出數據範圍，跳過。")
            continue

        # For demonstration, we'll analyze the segment from start_second to current_time_seconds
        # In a real application, you would implement proper ECG algorithms here.
        segment_start_sample = int(start_second * sample_rate)
        ecg_segment = ecg_data[segment_start_sample:current_sample_index]

        # Placeholder for actual ECG analysis
        vt_count = 0
        vt_duration_seconds = 0
        vf_count = 0
        vf_duration_seconds = 0
        vpc_count = 0
        q_wave_total_duration_seconds = 0
        heart_rate = np.nan # Placeholder for heart rate calculation

        if ecg_segment.size > 0:
            pass

        all_ecg_results.append({
            'time_point': timedelta(seconds=current_time_seconds),
            'vt_count': vt_count,
            'vt_duration_seconds': vt_duration_seconds,
            'vf_count': vf_count,
            'vf_duration_seconds': vf_duration_seconds,
            'vpc_count': vpc_count,
            'q_wave_total_duration_seconds': q_wave_total_duration_seconds,
            'heart_rate': heart_rate
        })
    return all_ecg_results

#-------------------------------------------------------------------------------------------------------------------------
# Main program entry
if __name__ == "__main__":
    acq_filepath = "/home/celeron0517/Desktop/實驗數據/20250430.acq"  # 替換為.acq檔案路徑

    # Read the file
    data_manager = read_acq_data(acq_filepath)

    if data_manager:
        # Get sampling rate
        sample_rate = data_manager.samples_per_second
        print(f"檔案採樣率: {sample_rate} Hz")

        # Get user input for the starting point
        while True:
            try:
                start_point_seconds = int(input("請輸入起始秒數 (例如: 0 代表檔案開頭): "))
                if start_point_seconds < 0:
                    print("起始秒數不能為負數。請重新輸入。")
                else:
                    break
            except ValueError:
                print("輸入無效。請輸入一個整數。")

        # Define time intervals relative to the start_point_seconds
        # First set of intervals: -5 min (before start_point), 0 (at start_point), +5, +10, +15, +20, +25, +30 min
        first_set_intervals_minutes = [-5, 0, 5, 10, 15, 20, 25, 30]
        first_set_intervals_seconds = [t * 60 for t in first_set_intervals_minutes]

        print("\n--- 第一組分析 (自起始點前5分鐘，及起始點與其後5, 10, 15, 20, 25, 30分鐘) ---")

        # Try to get blood pressure channel
        bp_channel = get_channel_data(data_manager, ['Blood Pressure', 'BP', 'Arterial Pressure'])
        if bp_channel:
            bp_results_first_set = calculate_blood_pressure_metrics(bp_channel, sample_rate, start_point_seconds, first_set_intervals_seconds)
            print("\n血壓分析結果 (第一組):")
            for res in bp_results_first_set:
                print(f"時間點: {res['time_point']}, 瞬時血壓最大值 (約0.5秒內): {res['instantaneous_bp_max']:.2f}, 瞬時血壓最小值 (約0.5秒內): {res['instantaneous_bp_min']:.2f}, 平均血壓 (自起始點到該點): {res['average_bp']:.2f}")
        else:
            print("未找到血壓通道，跳過血壓分析。")

        # Try to get ECG channel
        ecg_channel = get_channel_data(data_manager, ['ECG', 'Electrocardiogram', 'Lead II'])
        if ecg_channel:
            file = acq_filepath
            ecg_analysis_results_first_set = analyze_ecg_for_arrhythmias(file ,ecg_channel, sample_rate, start_point_seconds, first_set_intervals_seconds)
            print("\nECG心律不整分析結果 (第一組，待完成):")
            for res in ecg_analysis_results_first_set:
                print(f"時間點: {res['time_point']}, VT 次數: {res['vt_count']}, 總持續時間: {res['vt_duration_seconds']:.2f} 秒, VF 次數: {res['vf_count']}, 總持續時間: {res['vf_duration_seconds']:.2f} 秒, VPC 次數: {res['vpc_count']}, Q波總秒數: {res['q_wave_total_duration_seconds']:.2f} 秒, 心率: {res['heart_rate']:.2f} BPM")
        else:
            print("未找到ECG通道，跳過ECG分析。")

        # Get user input for the next starting point
        while True:
            try:
                next_start_point_seconds = int(input("\n請輸入下一個起始秒數 (例如: 1800 代表第30分鐘): "))
                if next_start_point_seconds < 0:
                    print("起始秒數不能為負數。請重新輸入。")
                else:
                    break
            except ValueError:
                print("輸入無效。請輸入一個整數。")

        # Second set of intervals: 0 (at next_start_point), +5, +10, +15, +30, +45, +60, +90, +120 min
        second_set_intervals_minutes = [0, 5, 10, 15, 30, 45, 60, 90, 120]
        second_set_intervals_seconds = [t * 60 for t in second_set_intervals_minutes]

        print(f"\n--- 第二組分析 (自 {next_start_point_seconds} 秒處，及該點後 5, 10, 15, 30, 45, 60, 90, 120 分鐘) ---")

        if bp_channel:
            bp_results_second_set = calculate_blood_pressure_metrics(bp_channel, sample_rate, next_start_point_seconds, second_set_intervals_seconds)
            print("\n血壓分析結果 (第二組):")
            for res in bp_results_second_set:
                print(f"時間點: {res['time_point']}, 瞬時血壓最大值 (約0.5秒內): {res['instantaneous_bp_max']:.2f}, 瞬時血壓最小值 (約0.5秒內): {res['instantaneous_bp_min']:.2f}, 平均血壓 (自起始點到該點): {res['average_bp']:.2f}")
        else:
            print("未找到血壓通道，跳過血壓分析。")

        if ecg_channel:
            ecg_analysis_results_second_set = analyze_ecg_for_arrhythmias(file, ecg_channel, sample_rate, next_start_point_seconds, second_set_intervals_seconds)
            print("\nECG心律不整分析結果 (第二組，待完成):")
            for res in ecg_analysis_results_second_set:
                print(f"時間點: {res['time_point']}, VT 次數: {res['vt_count']}, 總持續時間: {res['vt_duration_seconds']:.2f} 秒, VF 次數: {res['vf_count']}, 總持續時間: {res['vf_duration_seconds']:.2f} 秒, VPC 次數: {res['vpc_count']}, Q波總秒數: {res['q_wave_total_duration_seconds']:.2f} 秒, 心率: {res['heart_rate']:.2f} BPM")
        else:
            print("未找到ECG通道，跳過ECG分析。")