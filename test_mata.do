*！ description: do the local polynomial regression of Y on X and V with the grid of x and v
*! Meiting Wang
*! wangmeiting92@gmail.com
*! Oct 25, 2021


****** 参数设置
clear all
macro drop _all
cls
global n = 10*200
global h = 0.1 //bandwidth
global p = 3 //degree
global grid_num = 50
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


****** 调入 Mata 函数
do two_dimen_lpoly.mata
mata: mata mosave two_dimen_lpoly(), replace
lmbuild ltwo_dimen_lpoly, replace dir(.)


****** 使用 Mata 函数
timer clear 1
timer on 1
mata:
Y_var = st_data(.,"Y")
X_var = st_data(.,"X")
V_var = st_data(.,"V")
x_var = st_data(.,"x")
v_var = st_data(.,"v")
two_dimen_lpoly(Y_var,X_var,V_var,x_var,v_var,"$kernel",$h,$p)
end
timer off 1
timer list





