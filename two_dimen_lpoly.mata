mata:
real matrix function two_dimen_lpoly(
	real colvector Y_var, //dependent variable
	real colvector X_var, //independent variable 1
	real colvector V_var, //independent variable 2
	real colvector x_var, //grid point 1(the missing value must be at the end, if any.)
	real colvector v_var, //grid point 2(the missing value must be at the end, if any.)
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

	/* 得到没有缺漏值的 Y_raw X_raw V_raw x v，(Y_raw X_raw V_raw) 匹配，(x v)匹配 */
	Y_X_V_nomiss = (Y_var:!=.):&(X_var:!=.):&(V_var:!=.)
	x_v_nomiss = (x_var:!=.):&(v_var:!=.)
	Y_raw = select(Y_var,Y_X_V_nomiss)
	X_raw = select(X_var,Y_X_V_nomiss)
	V_raw = select(V_var,Y_X_V_nomiss)
	x = select(x_var,x_v_nomiss)
	v = select(v_var,x_v_nomiss)

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

	/* 得到符合 derivative 的 BETA(其列分别对应 beta00 beta10 beta01 ... betap0 .. beta0p) */
	multiplier = 1
	for (s=1; s<=p; s++) {
		multiplier = multiplier, (factorial(s..0):*factorial(0..s):/h^s)
	}
	BETA = BETA * diag(multiplier)

	/* 返回 BETA */
	return(BETA)
}
end
