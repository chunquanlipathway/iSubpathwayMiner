data<-read.table("C:/Users/熙/Documents/Tencent Files/934600892/FileRecv/data.xlsx",sep="\t",header=T)	#导入数据
re_t<-c()	#为T值构建向量
re_p<-c()	#为P值构建向量
for(i in 1:dim(data)[1])	#循环语句，i为变量，dim显示数据的行列数，dim（data)[1]表示到数据的最后一行，整条语句表示为第一行循环到最后一行
{
x<-as.numeric(data[i,2:791])	#对第2列达到第791列的值进行强制数字化处理，as.numeric为相应代码
y<-as.numeric(data[i,792:dim(data)[2]])	#对第792列达到最后一列的值进行强制数字化处理，as.numeric为相应代码，[2]代表列
re<-t.test(x,y)	#对x，y进行求值，对应的T值，P值
re_t<-c(re_t,re[[1]])	#求出相应的T值，将每次循环后的值放入相应的向量之中，[1]代表T值
re_p<-c(re_p,re[[3]])	#求出相应的P值，将每次循环后的值放入相应的向量之中，[3]代表P值
}
result<-cbind(as.character(data[,1]),re_t,re_p)	#输出结果，将T值P值放入矩阵当中，data[,1]表示数据所有行的一列，即基因名，as.character表示强制字符化
colnames(result)<-c("Gene's Names","T-Value","P-Value")	#更改矩阵列名，colnames为列名
write.table(result,"C:/Users/熙/Documents/Tencent Files/934600892/FileRecv/T值P值.txt",sep="\t",quote=F)	输出结果，quote=F表示吧结果中的“ ”去除