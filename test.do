****** 参数设置
clear all
macro drop _all
cls
global n = 10000
global h = 0.2 //bandwidth
global p = 3 //degree
global kernel "rectangle" //kernel function(the function written must be one of gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen, gaussian_m)
global seed = 123456


****** 生成测试数据
if "$seed" != "" {
	set seed $seed
}
set obs $n
gen X = rnormal()
gen V = rnormal()
gen Y = 1 + 5*X + 3*V + 2*X*V + rnormal()
range x -1 1 50
range v 0 1 50
order Y X V x v


****** lpoly2 命令的使用
lpoly2 Y X V, at(x v) kernel($kernel) bwidth($h) degree($p)
lpoly2 Y X V, at(x v) kernel($kernel) bwidth($h) degree($p) keep(beta00 beta10 beta01)
lpoly2 Y X V, at(x v) kernel($kernel) bwidth($h) degree($p) keep(beta00-beta01)
