

#Take in a fasta and an orftable and output a fastp
library(Biostrings)

translateFasta<-function(fi,fo,min.len=0){
	unlink(fo)
	fasta<-file(fi,'r')
	lastname<-""
	line<-readLines(fasta,1)
	while(length(line)>0){
		if(substr(line,1,1)=='>'){
			name<-strsplit(line," ")[[1]][1]
			name<-substr(name,2,nchar(name))
			if(nchar(lastname)>0){
				write.orfs(t.seq,lastname,fo)
				} 
			lastname<-name
			t.seq<-c()
			} else {
			t.seq<-c(t.seq,strsplit(tolower(line),"")[[1]])
			}
		line<-readLines(fasta,1)
		if(length(line)>0){if(line==""){line<-readLines(fasta,1)}}
		}
	write.orfs(t.seq,lastname,fo)
	close(fasta)
	}
	
write.orfs<-function(t.seq,name,fo){
	aa<-s.translate(t.seq)
	name<-paste(">",name,sep="")
	write(name,file=fo,append=TRUE)
	write(aa,file=fo,append=TRUE)
	}
	
s.translate<-function(t.seq){
	aa<-paste(sapply(seq(from=1,to=length(t.seq),by=3),function(q){
		cod<-t.seq[q:(q+2)]
		if('n' %in% cod | any(is.na(cod))){return('X')}
		return(as.character(GENETIC_CODE[[toupper(paste(cod,collapse=""))]]))
		}),collapse="")
	while(length(grep('XX',aa))>0){
		aa<-sub('XX','X',aa)
		}
	return(aa)
	}
