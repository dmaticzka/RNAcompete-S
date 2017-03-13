# pcws2_reg computes the piecewise linear regression -- with two pieces -- to (x,y), for any possible change point and chooses the one leading to the smallest least-square error.

pcws2_reg <- function(x, y)
{
	
	C <- length(x)
	ssBest = Inf
	for (c in 2:(C-1))
	{
		x1 <- x[1:c]
		y1 <- y[1:c]
		x2 <- x[c:C]
		y2 <- y[c:C]
		
		a1 <- sum((x1-mean(x1))*(y1-mean(y1)))/sum((x1-mean(x1))^2)
		b1 <- -a1 * mean(x1) + mean(y1)
		
		a2 <- sum((x2-mean(x2))*(y2-mean(y2)))/sum((x2-mean(x2))^2)
		b2 <- -a2 * mean(x2) + mean(y2)
		
		ss <- sum((a1*x1+b1-y1)^2) + sum((a2*x2+b2-y2)^2)
		
		if (ss < ssBest) 
		{
			ssBest <- ss
			cBest <- c
			a1Best <- a1
			a2Best <- a2
			b1Best <- b1
			b2Best <- b2
		}
	}
	
	return(list(c=cBest, a1=a1Best, b1=b1Best, a2=a2Best, b2=b2Best, residuals = c(a1*x1+b1-y1,a2*x2+b2-y2)))
}


# pcws3_reg computes the piecewise linear regression -- with three pieces -- to (x,y), for any possible change points and chooses the ones leading to the smallest least-square error.

pcws3_reg <- function(x, y)
{
	
	C <- length(x)
	ssBest = Inf
	for (c1 in 2:(C-2))
	{
		for (c2 in (c1+1):(C-1))
		{
			
			x1 <- x[1:c1]
			y1 <- y[1:c1]
			x2 <- x[c1:c2]
			y2 <- y[c1:c2] 
			x3 <- x[c2:C]
			y3 <- y[c2:C] 
			
			a1 <- sum((x1-mean(x1))*(y1-mean(y1)))/sum((x1-mean(x1))^2)
			b1 <- -a1 * mean(x1) + mean(y1)
			
			a2 <- sum((x2-mean(x2))*(y2-mean(y2)))/sum((x2-mean(x2))^2)
			b2 <- -a2 * mean(x2) + mean(y2)
			
			a3 <- sum((x3-mean(x3))*(y3-mean(y3)))/sum((x3-mean(x3))^2)
			b3 <- -a3 * mean(x3) + mean(y3)
			
			ss <- sum((a1*x1+b1-y1)^2) + sum((a2*x2+b2-y2)^2) + sum((a3*x3+b3-y3)^2)
			
			if (ss < ssBest) 
			{
				ssBest <- ss
				c1Best <- c1
				c2Best <- c2
				a1Best <- a1
				b1Best <- b1
				a2Best <- a2
				b2Best <- b2
				a3Best <- a3
				b3Best <- b3
			}
		}
	}
	return(list(c1=c1Best, c2=c2Best, a1=a1Best, b1=b1Best, a2=a2Best, b2=b2Best, a3=a3Best, b3=b3Best, residuals = c(a1*x1+b1-y1,a2*x2+b2-y2,a3*x3+b3-y3)))
}



xlog <- function(x) 
{
	xlog1d <- function (xi) if (xi == 0) 0 else (xi*log(xi))
	
	if (is.null(dim(x)))
	{
		return(sapply(x,xlog1d))
	}
	else
	{
		return(matrix(sapply(x,xlog1d),dim(x)))
	}
}

