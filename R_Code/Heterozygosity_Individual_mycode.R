count = 0
temp = NULL
for(i in 1:nrow(data))
{
	count = 0
	print(i)
	for(j in  seq(from=7,to=ncol(data),by=2))
	{
		if(as.character(data[i,j]) != as.character(data[i,j+1]))
		{
			count = count + 1;
		}
	}
	temp[i] = count
}
names(temp) = data[,2]