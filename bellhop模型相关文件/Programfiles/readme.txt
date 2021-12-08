cal_bellhop.m是计算Bellhop的主程序，包括写入环境文件、地形文件，计算声场，绘制声场等

Bellhop 计算海洋声场的结构&流程

流程：
输入环境文件、地形文件等→计算bellhop→Plot 声场

结构：

输入：
环境文件(.env) 
海面形状(.ati)	   默认平整海面
海底地形(.bty)	   需要根据实际地形输入
顶部反射系数(.trc) 可自己设置 也可以用默认
底部反射系数(.brc) 可自己设置 也可以用默认
声速剖面(.ssp)     可单独输入声速剖面文件，也可在环境文件中设置声速剖面

过程：
bellhop.exe   Fortran版本，已编译为windows可执行文件，可直接运行，根据输入文件计算声场 路径：atWin10_2020_11_4\windows-bin-20201102\bellhop.exe Fortran代码路径：atWin10_2020_11_4\Bellhop
bellhopM.m    matlab 版本  具体声线如何计算声场的物理计算方式可看这个 路径：atWin10_2020_11_4\Matlab\Bellhop\BellhopM.m

输出：
声场文件(.shd)
描述声线、本征声线文件(.ray)



Bellhop环境文件描述：
根据仿真要求依次设置
Title 		文件名
Frequency	声波频率
Nmedia 		介质层数（设1代表海底仅一层沉积层）
Options		海面边界条件设置  
	Option1 :声速剖面插值类型 C：三次多项式插值 N：N2型线性插值 S：三次样条插值 A：解析形式 
	Option2：顶部边界条件 V：真空 A：弹性半空间 R：绝对硬 F：由.trc文件输入 W：将界面反射系数写入.irc文件并作为底面边界条件读入 
	Option3: 衰减系数单位 N：Nepers/m F: dB/kmHz M:dB/m W:dB/wavelength Q:品质因数
Nmesh：进行数值计算的网格数，声学介质每波长10个点，弹性介质每波长20个点(计算bellhop默认为0) Sigma(m)：界面粗糙度，水平分层一般取0 Z(NSSP)：该层介质底面深度
SSP		声速剖面  （深度-纵波声速-横波声速-密度-纵波衰减系数-横波衰减系数）
Options         底部边界类型 V A R F 同顶部边界 *代表从.bty文件中读入地形文件
                底部深度 底部声速cb 横波声速 介质密度 纵波衰减 横波衰减
Number of src	声源个数
Source depth	声源深度，仅设置最浅和最深的深度，中间按照声源个数等间隔分布
Number of rcv depth	垂直方向接收个数 
Receiver depth	 （m）  垂直方向接收深度，仅设置最浅和最深的深度，中间按照接收个数等间隔分布   （这个决定了声场在竖直方向上的分布情况）
Number of rcv range	水平方向接收个数 
Receiver range	（km）  水平方向接收深度，仅设置最近和最袁的距离，中间按照接收个数等间隔分布   （这个决定了声场在水平方向上的分布情况）
Options		输出文件类型 Option1 :C = .shd声场 Option2 :G 几何波束
Number of beams	 计算声线数
Alpha 		take-off angle of src  声源出射角(deg)
Step zbox rbox   分别为计算步长(0代表自动计算步长，最深深度，最远距离)	

