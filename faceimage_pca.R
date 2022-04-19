library(rgdal)
training_images_folder = "C:/MS/FALL_2021/MULTIVARIATE/R_assignments/Faceimage_data/yalefaces_train/"
list_of_files = list.files(training_images_folder)

##################The value of n#################
n = length(list_of_files)
n
images_store = list()
for (i in 1:n) 
images_store[[i]] = readGDAL(paste(training_images_folder,list_of_files[i],sep=""))$band1
image_size_rows = 243
x11()
imagedata = matrix(images_store[[1]],nrow= image_size_rows,byrow=T)
image(t(imagedata[nrow(imagedata):1,]),col = grey(seq(0,1, length=256)))
############## value of P#######################
p = length(images_store[[1]]) 
p
x_data = matrix(0,n,p)
for (i in 1:n) x_data[i, ] = images_store[[i]]
x_bar = colMeans(x_data)
mu_hat = x_bar

################# plotting the avg face #############

x11()
imagedata = matrix(mu_hat,nrow= image_size_rows,byrow=T)
image(t(imagedata[nrow(imagedata):1,]),col = grey(seq(0,1,
length=256)),xlab="Average Face")
B = matrix(0,p,n)
for (i in 1:n)
{
 B[,i] = x_data[i, ] - x_bar
}
B_transpose_B = t(B)%*%B
eigen_results = eigen(B_transpose_B,only.values=F)
sum_eigs = sum(diag(B_transpose_B))
ELOV = array(0,n)
for (k in 1:n)
{
 ELOV[k] = 1 - sum(eigen_results$values[1:k])/sum_eigs
 print(paste(k," ",ELOV[k],sep=""))
}
################# pca using k = 8 #########################
k = 8
eigen_kvalues = eigen_results$values[1:k]
eigen_kvectors = eigen_results$vectors[,1:k]
eigvec_Sigma_hat_k = matrix(0,p,k)
for (i in 1:k) eigvec_Sigma_hat_k[ ,i] =
c(1/sqrt(eigen_kvalues[i]))*(B%*%matrix(eigen_kvectors[,i],ncol=1))
Ak_hat = t(eigvec_Sigma_hat_k)
x11()
imagedata = matrix(Ak_hat[3, ],nrow= image_size_rows,byrow=T)
image(t(imagedata[nrow(imagedata):1,]),col = grey(seq(0,1,length=256)),
xlab="Eigenface Number: (3)")
New_image = readGDAL("C:/MS/FALL_2021/MULTIVARIATE/R_assignments/Faceimage_data/yalefaces_test/subject01.happy")$band1
x11()
imagedata = matrix(New_image,nrow= image_size_rows,byrow=T)
image(t(imagedata[nrow(imagedata):1,]),col = grey(seq(0,1,length=256)), xlab=
"Test image")
PC_scores = Ak_hat%*%( matrix( New_image - mu_hat, ncol=1) )
Reconstructed_image = mu_hat + t(Ak_hat)%*%matrix(PC_scores,ncol=1)
x11()
imagedata = matrix(Reconstructed_image,nrow= image_size_rows,byrow=T)
image(t(imagedata[nrow(imagedata):1,]),col = grey(seq(0,1,length=256)),
xlab="Reconstructed test image with k = (8)
eigenfaces")


######################### pca using k = 38##########################

k = 38
eigen_kvalues = eigen_results$values[1:k]
eigen_kvectors = eigen_results$vectors[,1:k]
eigvec_Sigma_hat_k = matrix(0,p,k)
for (i in 1:k) eigvec_Sigma_hat_k[ ,i] =
c(1/sqrt(eigen_kvalues[i]))*(B%*%matrix(eigen_kvectors[,i],ncol=1))
Ak_hat = t(eigvec_Sigma_hat_k)
x11()
imagedata = matrix(Ak_hat[3, ],nrow= image_size_rows,byrow=T)
image(t(imagedata[nrow(imagedata):1,]),col = grey(seq(0,1,length=256)),
xlab="Eigenface Number: (3)")
New_image = readGDAL("C:/MS/FALL_2021/MULTIVARIATE/R_assignments/Faceimage_data/yalefaces_test/subject01.happy")$band1
x11()
imagedata = matrix(New_image,nrow= image_size_rows,byrow=T)
image(t(imagedata[nrow(imagedata):1,]),col = grey(seq(0,1,length=256)), xlab=
"Test image")
PC_scores = Ak_hat%*%( matrix( New_image - mu_hat, ncol=1) )
Reconstructed_image = mu_hat + t(Ak_hat)%*%matrix(PC_scores,ncol=1)
x11()
imagedata = matrix(Reconstructed_image,nrow= image_size_rows,byrow=T)
image(t(imagedata[nrow(imagedata):1,]),col = grey(seq(0,1,length=256)),
xlab="Reconstructed test image with k = (38)
eigenfaces")


