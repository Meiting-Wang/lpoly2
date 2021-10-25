****** 参数设置
clear all
macro drop _all
cls
global n = 10*200
global h = 0.1 //bandwidth
global p = 3 //degree
global grid_num = 50*200
global kernel "gaussian" //kernel function(the function written must be one of gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen, gaussian_m)
global seed = 123456


****** 生成测试数据
if "$seed" != "" {
	set seed $seed
}
set obs $n
gen X = rnormal()
gen V = rnormal()
gen Y = 1 + 5*X + 3*V + 2*X*V + rnormal()
range x -1 1 $grid_num
range v 0 1 $grid_num
order Y X V x v

****** lpoly2 命令的使用
timer clear 1
timer on 1
lpoly2 Y X V, at(x v) kernel($kernel) bwidth($h) degree($p)
timer off 1
timer list
