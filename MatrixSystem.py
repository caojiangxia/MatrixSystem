import numpy as np
import math
'''
输入： matrix为numpy格式
中间参数：
row,col为行数以及列数
L，R为待求矩阵
now_row,del_row为当前的主元行，和待消主元行
输出：
LU分解得到的
L矩阵
R矩阵
L*R结果
'''
def LU_Decomposition(matrix):
    row,col=matrix.shape
    if row!=col:#非方阵不可以进行LU分解
        print("we can't do LU decomposition! because the row is not equal the col")
        return
    L=np.zeros(shape=(row,col))
    for now_row in range(row):
        if matrix[now_row][now_row]==0:#当前行为非主元行，结束
            print("we can't do LU decomposition! because the row ",now_row," is not the main element")
            return
        L[now_row][now_row]=1
        for del_row in range(now_row+1,row):#使用主元行消去其他行
            L[del_row][now_row]=matrix[del_row][now_row]/matrix[now_row][now_row]#填充L矩阵的系数
            matrix[del_row]-=L[del_row][now_row]*matrix[now_row]#消去
    U=matrix
    print("L=")
    print(L)
    print("U=")
    print(U)
    print("L*U=")
    print(np.dot(L, U))
'''
输入： matrix为numpy格式
中间参数：
row,col为行数以及列数
Q，R为待求矩阵
now_col,del_col为当前的主元列，和已标准正交列
u表示当前列的引用
输出：
QR分解得到的
Q矩阵
R矩阵
Q*R结果
'''
def QR_Decomposition(matrix):
    row, col = matrix.shape
    R=np.zeros(shape=(col,col))
    matrix=matrix.T#将矩阵转置方便得到每一列
    for now_col in range(col):
        u=matrix[now_col]#取出当前列
        for del_col in range(now_col):
            R[del_col][now_col]=np.dot(u,matrix[del_col].T)#使用得到的标准基对当前列进行消去
            u-=R[del_col][now_col]*matrix[del_col]
        R[now_col][now_col]=math.sqrt(np.sum(np.square(u)))
        if math.sqrt(np.sum(np.square(u)))<0.00001:#若消去后为0向量，说明线性相关
            print("we can't do decomposition ,because the matrix linearly dependent")
            return
        u/=math.sqrt(np.sum(np.square(u)))
    Q=matrix.T
    print("Q=")
    print(Q)
    print("R=")
    print(R)
    print("Q*R=")
    print(np.dot(Q, R))
'''
输入： matrix为numpy格式
中间参数：
row,col为行数以及列数
Q，R为待求矩阵
now_col为当前的主元列
u表示当前列的引用
r表示生成的变换矩阵
输出：
QR分解得到的
Q矩阵
R矩阵
Q*R结果
'''
def OR_Decomposition(matrix):
    row, col = matrix.shape
    Q=np.eye(row)
    for now_col in range(col-1):
        u=np.array([matrix.T[now_col:,now_col:][0]])#提取当前列
        u[0][0]=u[0][0]-math.sqrt(np.sum(np.square(u)))
        r = np.eye(row)
        r[now_col:, now_col:]=r[now_col:,now_col:]-2*np.dot(u.T,u)/np.dot(u,u.T)[0][0]#计算r矩阵
        matrix=np.dot(r,matrix)#对当前列标准化
        if math.fabs(matrix[now_col+1][now_col+1])<0.00001:#若为0向量，说明线性相关
            print("we can't do decomposition ,because the matrix linearly dependent")
            return
        Q=np.dot(r,Q)
    Q=Q.T
    R=matrix
    print("Q=")
    print(Q)
    print("R=")
    print(R)
    print("Q*R=")
    print(np.dot(Q, R))
'''
输入： matrix为numpy格式
中间参数：
row,col为行数以及列数
Q，R为待求矩阵
now_row为当前的主元列
res旋转时的分母
r表示生成的旋转矩阵
输出：
QR分解得到的
Q矩阵
R矩阵
Q*R结果
'''
def GR_Decomposition(matrix):
    row, col = matrix.shape
    Q = np.eye(row)
    for now_row in range(col-1):#提取当前列
        res=matrix[now_row][now_row]*matrix[now_row][now_row]
        for del_row in range(now_row+1,row):
            r = np.eye(row)
            res+=matrix[del_row][now_row]*matrix[del_row][now_row]#分母部分
            if res < 0.00001:#若为0向量，无法构造sin和cos
                print("we can't do decomposition ,because the matrix linearly dependent")
                return
            r[now_row][now_row]=0
            r[del_row][del_row] = 0#进行旋转
            r[now_row][now_row] = matrix[now_row][now_row]/math.sqrt(res)
            r[del_row][del_row] = matrix[now_row][now_row]/math.sqrt(res)
            r[now_row][del_row] = matrix[del_row][now_row]/math.sqrt(res)
            r[del_row][now_row] = -matrix[del_row][now_row]/math.sqrt(res)
            matrix = np.dot(r, matrix)
            Q = np.dot(r, Q)
        if math.fabs(matrix[now_row+1][now_row+1]) < 0.00001:#若为0向量，说明线性相关
            print("we can't do decomposition ,because the matrix linearly dependent")
            return
    Q = Q.T
    R = matrix
    print("Q=")
    print(Q)
    print("R=")
    print(R)
    print("Q*R=")
    print(np.dot(Q, R))
def Decomposition(operate,matrix):
    matrix=np.array(matrix,dtype=np.float64)
    row,col=matrix.shape
    if col > row :
        print("col > row ,the martix can't do ",operate," decomposition!")
        return
    if operate=="LU":
        print("Begin LU decomposition:")
        LU_Decomposition(matrix)
    if operate=="QR":
        print("Begin QR(Gram-Schmidt) decomposition:")
        QR_Decomposition(matrix)
    if operate=="OR":
        print("Begin OR(Householder reduction) decomposition:")
        OR_Decomposition(matrix)
    if operate=="GR":
        print("Begin GR(Givens reduction) decomposition:")
        GR_Decomposition(matrix)
if __name__ == '__main__':
    while 1:
        n=int(input())#输入行数
        if n==-1:
            break
        matrix=[]
        for i in range(n):#输入矩阵
            res=input()
            res=res.strip().split(" ")
            for j in range(len(res)):
                res[j]=int(res[j])
            matrix.append(res)
        operate=input()#输入操作
        Decomposition(operate,matrix)#开始分解
'''
example:

LU:
3
2 2 2
4 7 7
6 18 22

answer :
L:
1 0 0
2 1 0
3 4 1
U:
2 2 2
0 3 3
0 0 4

QR:
3
0 -20 -14
3 27 -4
4 11 -2
Q:
1/25
0 -20 -15
15 12 -16
20 -9 12
R:
5 25 -4
0 25 10
0 0 10


QR
4
4 -3 4
2 -14 3
-2 14 0
1 -7 15


QR
3
4 2 -2 1
-3 -14 14 7
4 3 0 15
'''