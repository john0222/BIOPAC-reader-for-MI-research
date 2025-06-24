
import numpy as np
import scipy.signal as signal
from collections import deque
import matplotlib.pyplot as plt
import bioread # 導入 bioread 函式庫

from pylab import mpl
mpl.rcParams["font.sans-serif"] = ["SimHei"]
mpl.rcParams["axes.unicode_minus"] = False
# RatECGAnalyzer class
class RatECGAnalyzer:
    """
    rr interval analyzer，Pan-Tompkins algo for rat heart
    """
    def __init__(self, sfreq, highpass_cutoff=5.0, lowpass_cutoff=40.0,
                 integration_window_ms=100, threshold_factor=0.25, noise_threshold_factor=0.5,
                 refractory_period_ms=150):
        """
        initalization on the analyzer

            sfreq (float): 採樣頻率 (Hz)。
            highpass_cutoff (float): 帶通濾波器的高通截止頻率 (Hz)。
                                     大鼠心率較快，高通可設為5-10 Hz。
            lowpass_cutoff (float): 帶通濾波器的低通截止頻率 (Hz)。
                                    大鼠心率較快，低通可設為40-60 Hz。
            integration_window_ms (int): 移動平均濾波器窗口大小 (毫秒)。
                                         應與R波寬度相近，大鼠可設定為80-120 ms。
            threshold_factor (float): R波檢測閾值的因子。
            noise_threshold_factor (float): 雜訊閾值的因子。
            refractory_period_ms (int): 不應重複檢測R波的不應期 (毫秒)。
                                        大鼠心率快，不應期可設為100-200 ms。
        """
        self.sfreq = sfreq
        self.nyquist = 0.5 * sfreq
        self.highpass_cutoff = highpass_cutoff
        self.lowpass_cutoff = lowpass_cutoff
        self.integration_window = int(integration_window_ms * sfreq / 1000)
        self.threshold_factor = threshold_factor
        self.noise_threshold_factor = noise_threshold_factor
        self.refractory_period_samples = int(refractory_period_ms * sfreq / 1000)

        # init of algo status
        self.peak_values = deque(maxlen=int(sfreq * 2))  # r prak nearest
        self.noise_values = deque(maxlen=int(sfreq * 2)) # noise peak nearest
        self.recent_peak_location = -self.refractory_period_samples # r interval location
        # thereshold initalization
        self.r_peak_threshold = 0.0
        self.noise_peak_threshold = 0.0

    def _bandpass_filter(self, data):
        
        low = self.highpass_cutoff / self.nyquist
        high = self.lowpass_cutoff / self.nyquist
        b, a = signal.butter(1, [low, high], btype='band')
        return signal.filtfilt(b, a, data)

    def _differentiator(self, data):
     

        # 使用[1, 2, 0, -2, -1] / 8 的近似，考慮到大鼠心電圖頻率
        diff_filter = np.array([1, 2, 0, -2, -1]) / 8.0 * self.sfreq
        return signal.convolve(data, diff_filter, mode='same')

    def _squaring(self, data):
 
        return data**2

    def _moving_window_integration(self, data):
  
        return np.convolve(data, np.ones(self.integration_window) / self.integration_window, mode='same')

    def _update_thresholds(self, current_peak_val, is_r_peak=True):

        if is_r_peak:
            self.peak_values.append(current_peak_val)
        else:
            self.noise_values.append(current_peak_val)

        # slide window refresh 
        if len(self.peak_values) > 0:
            self.r_peak_threshold = np.mean(self.peak_values) * self.threshold_factor
        if len(self.noise_values) > 0:
            self.noise_peak_threshold = np.mean(self.noise_values) * self.noise_threshold_factor

        # make sure that r wave thereshold higher then noise 
        if self.r_peak_threshold < self.noise_peak_threshold * 1.2: # mild raise of ratio on r / noise 
            self.r_peak_threshold = self.noise_peak_threshold * 1.2

    def analyze_ecg(self, ecg_signal):
        """
        分析單通道ECG訊號，檢測R波並計算R-R間期。

        參數:
            ecg_signal (np.array): 一維的ECG訊號數據。

        返回:
            tuple: 包含以下元素的元組：
                - r_peak_indices (list): 檢測到的R波在原始訊號中的索引。
                - rr_intervals_ms (list): R-R間期列表 (毫秒)。
        """
        if not isinstance(ecg_signal, np.ndarray):
            ecg_signal = np.array(ecg_signal)

        # Pan-Tompkins algo 
        # bandpass filter 
        filtered_signal = self._bandpass_filter(ecg_signal)

        # diff
        differentiated_signal = self._differentiator(filtered_signal)

        # square 
        squared_signal = self._squaring(differentiated_signal)

        # slide windows 
        integrated_signal = self._moving_window_integration(squared_signal)

        r_peak_indices = []
        rr_intervals_ms = []

        peak_candidates = [] # to store the potential peak 
        search_back_buffer = [] # search back 

        for i in range(1, len(integrated_signal)):
            # raise dectect 
            if integrated_signal[i] > integrated_signal[i-1]:
                # descent dectect (potential peak)
                if i + 1 < len(integrated_signal) and integrated_signal[i] > integrated_signal[i+1]:
                    current_peak_val = integrated_signal[i]
                    current_peak_idx = i

                    # no reaction period check 
                    if (current_peak_idx - self.recent_peak_location) < self.refractory_period_samples:
                        # if still in period, ignore
                        search_back_buffer = [] # flush search back buffer
                        continue

                    # higher then current r wave 
                    if current_peak_val >= self.r_peak_threshold:
                        r_peak_indices.append(current_peak_idx)
                        if len(r_peak_indices) > 1:
                            rr_intervals_ms.append((r_peak_indices[-1] - r_peak_indices[-2]) * 1000 / self.sfreq)
                        self._update_thresholds(current_peak_val, is_r_peak=True)
                        self.recent_peak_location = current_peak_idx
                        search_back_buffer = [] # r detected,flush buffer 

                    # noise peak is higher than noise peak stored but lower than r peak stored 
                    elif current_peak_val >= self.noise_peak_threshold:
                        search_back_buffer.append((current_peak_idx, current_peak_val))
                        self._update_thresholds(current_peak_val, is_r_peak=False) # refresh noise as new thereshold 

                    # current peak lower than noise peak 
                    else:
                        search_back_buffer = [] # flush buffer
                        self._update_thresholds(current_peak_val, is_r_peak=False) # refresh noise as new thereshold 

            # 搜尋回溯機制：如果長時間沒有檢測到R波，回溯檢查緩衝區中的潛在峰值
            # 簡單判斷：如果超過預設時間沒有R波，且緩衝區中有峰值
            # 實際應用中可能需要更複雜的判斷，例如判斷平均心率等

            if (i - self.recent_peak_location) > self.sfreq * 0.8: # 約0.8s no r wave detected (基於大鼠心率，假設最低120 BPM)
                if len(search_back_buffer) > 0:
                    # search highest peak in buffer 
                    best_candidate_idx, best_candidate_val = max(search_back_buffer, key=lambda item: item[1])

                    # 僅當最高峰值顯著高於雜訊閾值時才考慮為R波
                    if best_candidate_val >= self.noise_peak_threshold * 1.5: # 提高回溯時的閾值，降低誤判
                        r_peak_indices.append(best_candidate_idx)
                        if len(r_peak_indices) > 1:
                            rr_intervals_ms.append((r_peak_indices[-1] - r_peak_indices[-2]) * 1000 / self.sfreq)
                        self._update_thresholds(best_candidate_val, is_r_peak=True)
                        self.recent_peak_location = best_candidate_idx
                        search_back_buffer = [] # r deccted, flush buffer

        return r_peak_indices, rr_intervals_ms
    

#-----------------------------------------------------------------------------------------
    def main(file_path):
        # FILE PATH
        acq_file_path = file_path
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
                    plt.title(f'大鼠ECG R波檢測結果 (通道: {ecg_channel.name})')
                    plt.xlabel('時間 (秒)')
                    plt.ylabel('振幅')
                    plt.grid(True)
                    plt.legend()
                    plt.tight_layout()
                    plt.show()

                    # RR interval plot
                    if rr_intervals_ms:
                        plt.figure(figsize=(10, 5))
                        plt.hist(rr_intervals_ms, bins=20, edgecolor='black')
                        plt.title('R-R間期分佈')
                        plt.xlabel('R-R間期 (毫秒)')
                        plt.ylabel('頻率')
                        plt.grid(True)
                        plt.tight_layout()
                        plt.show()

                        # calculate average bpm
                        average_rr_ms = np.mean(rr_intervals_ms)
                        average_hr_bpm = 60000 / average_rr_ms if average_rr_ms > 0 else 0
                        print(f"平均R-R間期: {average_rr_ms:.2f} 毫秒")
                        print(f"平均心率: {average_hr_bpm:.2f} BPM")
                    else:
                        print("沒有足夠的R-R間期數據來繪製直方圖或計算平均心率。")

                else:
                    print("在提供的 .acq 檔案中未找到任何可用的數據通道。")
            else:
                print(f"錯誤: 檔案 '{acq_file_path}' 中沒有任何通道數據。")

        except FileNotFoundError:
            print(f"錯誤: 找不到檔案 '{acq_file_path}'。請檢查路徑是否正確。")
        except Exception as e:
            print(f"讀取或處理檔案時發生錯誤: {e}")
        return r_peak_indices,rr_intervals_ms
