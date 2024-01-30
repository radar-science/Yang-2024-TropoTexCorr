Instruction:
   1) 主函数名称：大写字母开头；子函数名称：小写字母开头
   2) TextureCorrelation.m：仅包含基于纹理相关的对流层延迟矫正方法；TextureCorrelation_CompareVersion.m: 包含对比方法的代码
   3) 可调参数：对应代码行后用**注释。有的可调参数在子函数中，包含可调参数的子函数在主函数被调用行也以**注释。子函数中的参数一般情况下不用改变，子函数中的可调参数要么在程序前几行集中列出，要么同样以**注释。
Input：
   1) timeseries_ramp_demErr.h5 (must): 去除相位斜坡和DEM误差相位后相位时间序列，通过mintpy生成
   2) height.tif (must): 高程数据，通过save_gdal.py geo_geometryRadar.h5 -d height --of GTiff代码生成，其中geo_geometryRadar.h5是mintpy文件
   3) temporalCoherence.h5 (must): 时间相干性，用于生成掩膜，通过mintpy生成
   4) velocity.h5 (optional): 变速率数据，用于生成掩膜，通过mintpy生成
   5) timeseries_ERA5_ramp_demErr.h5 (optional): 去除相位斜坡和DEM误差相位并用ERA5矫正对流层延迟后的相位时间序列，用于对比和参考，通过mintpy生成
Output (all in .mat format)：
   1) phase_ts_deramp.m (compare)：Original phase time series (serve as input and comparison)
   2) phase_ts_atm.m (compare): Original phase time series after ERA5 correction.
   3) phase_ts_LLF.m (compare)：Original phase time series after local linear fitting method correction.
   4) phase_ts_GLF.m (compare)：Original phase time series after global linear fitting method correction.
   5) phase_ts_HTC_low.m: Original phase time series after texture correlation method (low-resolution module) correction.
   6) phase_ts_HTC_high.m: Original phase time series after texture correlation method (high-resolution module) correction.
   7) phase_ts_HTC.m：phase_ts_HTC.m can be assigned as phase_ts_HTC_low.m or phase_ts_HTC_high.m.
