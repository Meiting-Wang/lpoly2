* Description: Do two-dimensional local polynomial regression.
* Author: Meiting Wang, Ph.D. Candidate, Institute for Economic and Social Research, Jinan University
* Email: wangmeiting92@gmail.com
* Created on Oct 23, 2021


*-------------------------主程序-----------------------------
program define lpoly2
version 16.0

syntax varlist(min=3 max=3 numeric) , at(varlist numeric min=2 max=2) [KERnel(string) BWidth(numlist max=1 >0) Degree(numlist max=1 >=0 integer) keep(string)]

* 默认值(设定 degree keep kernel bwidth 的默认值)
if "`kernel'" == "" {
	local kernel "gaussian"
}
if "`bwidth'" == "" {
	local bwidth = 0.1
}
if "`degree'" == "" {
	local degree = 3
}
if "`keep'" == "" {
	local keep "beta00"
	if `degree' >= 1 {
		forvalues s = 1/`degree' {
			forvalues a = `s'(-1)0 {
				local b = `s'-`a'
				local keep "`keep' beta`a'`b'"
			}
		}
	}
}

* 错误输入识别
if ~ustrregexm("`keep'","(^beta\d{2}(\s+beta\d{2})*$)|(^beta\d{2}-beta\d{2}$)") {
	dis as error "Option {opt keep} synatx error"
	error 9999
} //保证 keep 输入的是 "beta00 beta10 ..." 或 "beta00-beta01" 的类似语句

* 获得因变量, 自变量和格点变量
tokenize `varlist'
local Y "`1'"
local X "`2'"
local V "`3'"
tokenize `at'
local x "`1'"
local v "`2'"

if ("`x'"=="`X'")|("`x'"=="`V'")|("`v'"=="`X'")|("`v'"=="`V'") {
	dis as error "Variables in {opt at(varlist)} and {it:xvars} can not be the same."
	error 9999
} //以保证 at(varlist) 和 xvars 中的变量不一样
if "`x'" == "`v'" {
	dis as error "Variables in {opt at(varlist)} can not be the same"
	error 9999
} //保证填入 at(varlist) 的变量不一样


* 使用子程序中的 two_dimen_lpoly() 函数，计算 two dimensional 的 local polynomial regression, 并保存所设定的变量
qui ds
local varlist_old "`r(varlist)'"
mata: two_dimen_lpoly("`Y'","`X'","`V'","`x'","`v'","`kernel'",`bwidth',`degree')
keep `varlist_old' `keep'


* 输入参数展示
dis as text "Parameters:"
dis _col(6) as text "bwidth = {result:`bwidth'}"
dis _col(6) as text "degree = {result:`degree'}"
dis _col(6) as text "kernel = {result:`kernel'}"
dis _col(6) _s(2) as text "keep = {result:`keep'}"

end



*--------------------------子程序-----------------------------
/* 能计算 two dimensional 的 local polynomial regression，并把结果储存在 Stata 数据集的对应变量中(beta00 beta10 beta01 ... betap0 .. beta0p) */
version 16.0
mata:
void function two_dimen_lpoly(
	string scalar Y_var, //dependent variable
	string scalar X_var, //independent variable 1
	string scalar V_var, //independent variable 2
	string scalar x_var, //grid point 1(the missing value must be at the end, if any.)
	string scalar v_var, //grid point 2(the missing value must be at the end, if any.)
	string scalar kernel, //kernel function(the function written must be one of the series of functions specified in the program)
	real scalar h, //bandwidth(bandwidth needs to be greater than 0)
	real scalar p) //degree(degree needs to be a non-negative integer)
{
	/* 错误信息判断 */
	if (h<=0) {
		printf("{error:Bandwidth needs to be greater than 0}\n")
		exit(9999)
	} //保证 h>0
	if ( (p<0) | (mod(p,1)) ) {
		printf("{error:Degree should be a non-negative integer}\n")
		exit(9999)
	} //保证 p 为非负整数

	/* 得到没有缺漏值的 Y_raw X_raw V_raw x v，且 listwise(Y_raw X_raw V_raw) 和 listwise(x v) */
	Y_raw = st_data(.,Y_var)
	X_raw = st_data(.,X_var)
	V_raw = st_data(.,V_var)
	x = st_data(.,x_var)
	v = st_data(.,v_var)
	Y_X_V_nomiss = (Y_raw:!=.):&(X_raw:!=.):&(V_raw:!=.)
	x_v_nomiss = (x:!=.):&(v:!=.)
	Y_raw = select(Y_raw,Y_X_V_nomiss)
	X_raw = select(X_raw,Y_X_V_nomiss)
	V_raw = select(V_raw,Y_X_V_nomiss)
	x = select(x,x_v_nomiss)
	v = select(v,x_v_nomiss)

	/* 得到 q */
	q = rows(x)

	/* 计算BETA */
	BETA = J(q,(p+2)*(p+1)/2,.) //设置结果储存矩阵
	for (i=1; i<=q; i++) {
		X_std_raw = (X_raw:-x[i,1])/h
		V_std_raw = (V_raw:-v[i,1])/h

		/* 依据不同的核函数对 k_index X_std V_std k 进行不同的处理(核函数的类别有: gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen, gaussian_m) */
		if (kernel=="gaussian") {
			k_index = J(rows(X_std_raw),1,1)
			X_std = select(X_std_raw,k_index)
			V_std = select(V_std_raw,k_index)
			k = normalden(X_std):*normalden(V_std)
		}
		else if (kernel=="epanechnikov") {
			k_index = (abs(X_std_raw):<sqrt(5)):&(abs(V_std_raw):<sqrt(5))
			X_std = select(X_std_raw,k_index)
			V_std = select(V_std_raw,k_index)
			k = (3/4*(1:-X_std:^2/5)/sqrt(5)):*(3/4*(1:-V_std:^2/5)/sqrt(5))
		}
		else if (kernel=="epan2") {
			k_index = (abs(X_std_raw):<1):&(abs(V_std_raw):<1)
			X_std = select(X_std_raw,k_index)
			V_std = select(V_std_raw,k_index)
			k = (3/4*(1:-X_std:^2)):*(3/4*(1:-V_std:^2))
		}
		else if (kernel=="biweight") {
			k_index = (abs(X_std_raw):<1):&(abs(V_std_raw):<1)
			X_std = select(X_std_raw,k_index)
			V_std = select(V_std_raw,k_index)
			k = (15/16*(1:-X_std:^2):^2):*(15/16*(1:-V_std:^2):^2)
		}
		else if (kernel=="cosine") {
			k_index = (abs(X_std_raw):<0.5):&(abs(V_std_raw):<0.5)
			X_std = select(X_std_raw,k_index)
			V_std = select(V_std_raw,k_index)
			k = (1:+cos(2*pi()*X_std)):*(1:+cos(2*pi()*V_std))
		}
		else if (kernel=="rectangle") {
			k_index = (abs(X_std_raw):<1):&(abs(V_std_raw):<1)
			X_std = select(X_std_raw,k_index)
			V_std = select(V_std_raw,k_index)
			k = (J(rows(X_std),1,0.5)):*(J(rows(V_std),1,0.5))
		}
		else if (kernel=="triangle") {
			k_index = (abs(X_std_raw):<1):&(abs(V_std_raw):<1)
			X_std = select(X_std_raw,k_index)
			V_std = select(V_std_raw,k_index)
			k = (1:-abs(X_std)):*(1:-abs(V_std))
		}
		else if (kernel=="parzen") {
			k_index = (abs(X_std_raw):<=1):&(abs(V_std_raw):<=1)
			X_std = select(X_std_raw,k_index)
			V_std = select(V_std_raw,k_index)

			k1 = J(rows(X_std),1,.)
			cond1 = selectindex(abs(X_std):<=0.5)
			cond2 = selectindex(abs(X_std):>0.5)
			k1[cond1] = 4/3:-8*X_std[cond1]:^2+8*abs(X_std[cond1]):^3
			k1[cond2] = 8*(1:-abs(X_std[cond2])):^3/3

			k2 = J(rows(V_std),1,.)
			cond1 = selectindex(abs(V_std):<=0.5)
			cond2 = selectindex(abs(V_std):>0.5)
			k2[cond1] = 4/3:-8*V_std[cond1]:^2+8*abs(V_std[cond1]):^3
			k2[cond2] = 8*(1:-abs(V_std[cond2])):^3/3

			k = k1:*k2
		}
		else if (kernel=="gaussian_m") {
			k_index = (abs(X_std_raw):<=6):&(abs(V_std_raw):<=6)
			X_std = select(X_std_raw,k_index)
			V_std = select(V_std_raw,k_index)

			k1 = J(rows(X_std),1,.)
			cond1 = selectindex(abs(X_std):<=5)
			cond2 = selectindex(abs(X_std):>5)
			k1[cond1] = normalden(X_std[cond1])
			k1[cond2] = normalden(5)*(4*(6:-abs(X_std[cond2])):^5-6*(6:-abs(X_std[cond2])):^4+3*(6:-abs(X_std[cond2])):^3)

			k2 = J(rows(V_std),1,.)
			cond1 = selectindex(abs(V_std):<=5)
			cond2 = selectindex(abs(V_std):>5)
			k2[cond1] = normalden(V_std[cond1])
			k2[cond2] = normalden(5)*(4*(6:-abs(V_std[cond2])):^5-6*(6:-abs(V_std[cond2])):^4+3*(6:-abs(V_std[cond2])):^3)

			k = k1:*k2
		}
		else {
			printf("{error:Wrong kernel function}\n")
			exit(9999)
		}

		/* 获得 Y Z */
		Y = select(Y_raw,k_index)
		Z = J(rows(X_std),1,1)
		for (s=1; s<=p; s++) {
			Z = Z,(mm_expand(X_std,1,s+1):^(s..0)):*(mm_expand(V_std,1,s+1):^(0..s))
		}

		/* 计算单组 x v 下的beta，并用其填充 BETA 矩阵的第 i 行 */
		beta = (cholinv(cross(Z,k,Z)))*(cross(Z,k,Y)) //cross(Z,k,Y) = Z'diag(w)Y
		BETA[i,.] = beta'
	}

	/* 得到符合 derivative 的 BETA */
	multiplier = 1
	for (s=1; s<=p; s++) {
		multiplier = multiplier, (factorial(s..0):*factorial(0..s):/h^s)
	}
	BETA = BETA * diag(multiplier)

	/* 将 BETA 转化为 Stata 内存中的变量(beta00 beta10 beta01 ... betap0 .. beta0p) */
	varnames = "beta00"
	for (s=1; s<=p; s++) {
		varnames = varnames,(J(1,s+1,"beta"):+strofreal(s..0):+strofreal(0..s))
	}
	st_store((1,rows(BETA)), st_addvar("double",varnames), BETA)
}
end
