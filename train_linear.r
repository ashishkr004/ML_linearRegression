#!/usr/bin/env Rscript
# your code must have this first line.
# Test code for linear regression part goes here

#function to get transpose of matrix
matrixTranspose<-function(matrixM,print="True"){
	elementOfTrans<-c()
	for (i in 1:nrow(matrixM)){
		for (j in 1:ncol(matrixM)){
			elementOfTrans<-append(elementOfTrans,matrixM[i,j])
		}
	}
	transposeOfMatrixX<-matrix(elementOfTrans,nrow=ncol(matrixM),ncol=nrow(matrixM))
	return(transposeOfMatrixX)
}

#matrix multiplication function
matrixMultiplication<-function(matrixA,matrixB,print="True"){
	if(ncol(matrixA)==nrow(matrixB)){
		elementOfAB<-c()
		for (col in 1:ncol(matrixB)){
			for (row in 1:nrow(matrixA)){
				num<-0
				for (i in 1:ncol(matrixA)){
         			num<-num+(matrixA[row,i]*matrixB[i,col])
				}
				elementOfAB<-append(elementOfAB,num)
			}
		}
		matrixAB<-matrix(elementOfAB,nrow=nrow(matrixA),ncol=ncol(matrixB))
		return(matrixAB)
	}
	else{
		return(0)
	}
}

#matrix subtraction
matrixSubtraction<-function(matrixA,matrixB,print="True"){
	if(nrow(matrixA)==nrow(matrixB)&ncol(matrixA)==ncol(matrixB)){
		for (i in 1:nrow(matrixA)){
			for (j in 1:ncol(matrixA)){
				matrixA[i,j]<-matrixA[i,j]-matrixB[i,j]
			}
		}
		return(matrixA)
	}
	else{
		return(0)
	}
}

#matrix inverse
matrixInverse<-function(matrixM,print="True"){
	if(nrow(matrixM)==ncol(matrixM)){
		return(solve(matrixM))
	}
	else{
		return(0)
	}
}


# Row echelon form
reducedRowEchelonForm<-function(matrix,print="True"){
    stopifnot(is.numeric(matrix))
    if (!is.matrix(matrix))
        stop("Input parameter must be a matrix.")

    numRows <- nrow(matrix)
    numCols <- ncol(matrix)-1

    dPal <- 0.00001
	j <- 1

    for (i in 1:numCols){     # for each columns of matrix.
    	pivot <- matrix[j:numRows,i]
    	pivot <- abs(pivot)
    	pivot <- which.max(pivot)
        pivot <- j+pivot-1

        m <- abs(matrix[pivot, i])

        if (m <= dPal) {
            matrix[j:numRows, i] <- 0  # zeros(nr-r+1, 1)
        }
        else {
            matrix[c(pivot, j), i:numCols] <- matrix[c(j, pivot), i:numCols]
            matrix[j, i:numCols] <- matrix[j, i:numCols] / matrix[j, i]
            if (j == 1){
                dx <- c((j+1):numRows)
            }
            else if (j == numRows){
                dx <- c(1:(j-1))
            }
            else {
                dx <- c(1:(j-1), (j+1):numRows)
            }
            matrix[dx, i:numCols] <- matrix[dx, i:numCols] - matrix[dx, i, drop=FALSE] %*% matrix[j, i:numCols, drop=FALSE]
            if(j == numRows){
            	break
            }
            j <- j+1
        }
    }
    matrix[abs(matrix) < dPal] <- 0
    vectorC <- c(1:numCols)
    matrix <- matrix[,vectorC]
    return(matrix)
}

# Column span matrix N of matrix M
spanMatrix<-function(matrixM,print="True"){
	rrefOfMatrixM<-reducedRowEchelonForm(matrixM)		# get reduced row echolon form of matrix M
	vecOfCols<-c()		#list of elements of columns which are independent
	
	#for each row of reduce echolon form matrix check each column is value is not equal to zero. 
	for(i in 1:nrow(rrefOfMatrixM)){
		for(j in 1:ncol(rrefOfMatrixM)){
			if(rrefOfMatrixM[i,j]==1){
				vecOfCols<-c(vecOfCols,j)		#column j is independent
				break
			}
			else if(rrefOfMatrixM[i,j]!=0){
				break
			}
		}
	}

	matrixN<-matrixM[,vecOfCols]		# Construct the new matrix N from all independents columns of matrix M.
	return(matrixN)
}

#pseudo inverse matrix function -> B=inverse of(tranpose(X)*X) * X.
pseudoInverseMatrix<-function(matrixM,print="True"){
	transposeOfMatrixM<-matrixTranspose(matrixM)
 	return(matrixMultiplication(matrixInverse(matrixMultiplication(transposeOfMatrixM,matrixM)),transposeOfMatrixM))
 }

# Gaussian elimination
eliminationFeature<-function(matrixM,print="True"){
	tempSpanMatrix<-spanMatrix(matrixM)
	return(pseudoInverseMatrix(tempSpanMatrix))
}



#main function
args<-commandArgs(TRUE)
data<-read.csv("input1.csv",header=FALSE, sep=",")		#read csv file and separate a data by comma.
matrix_A<-as.matrix(data)		#convert a data in matrix form.

# Construction of matrix X
listOfCols<-seq(1, ncol(matrix_A)-1, 1)
matrix_X<-matrix_A[,listOfCols]
columnOne<-matrix(1,nrow=nrow(matrix_X),ncol=1)
matrix_X<-cbind(columnOne,matrix_X)

# Construction of matrix Y
matrix_Y<-matrix_A[,ncol(matrix_A)]
matrix_Y<-matrix(matrix_Y,nrow=nrow(matrix_A),ncol=1)

#beta
beta<-matrixMultiplication(pseudoInverseMatrix(matrix_X),matrix_Y)

write.table(beta, file = "modelfile.csv", append = FALSE, quote = FALSE, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE, qmethod = c("escape", "double"),fileEncoding = "")












