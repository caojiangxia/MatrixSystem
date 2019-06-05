# 国科大-李保滨-矩阵大作业

## 程序说明

本程序基于python3.6,代码在MatrixSystem.py文件中,辅助使用numpy和math包，numpy中仅仅使用了二维矩阵的乘法、转置、内积。math包中仅仅使用sqrt方法。实现了LU分解、QR分解的三种方法("QR","OR","GR")

LU表示使用guass消元方式进行分解

QR表示使用Gram-Schmidt方式进行分解

OR表示使用Householder reduction方式进行分解

GR表示使用Givens reduction方式进行分解

### 运行说明

输入：
第一行输入n，表示矩阵行数。当n为-1时程序结束

第1行至第n+1行输入矩阵每行元素，每行m个数。

第n+2行输入分解方式("LU","QR","OR","GR"其中一个)
