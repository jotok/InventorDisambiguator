
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% This file rearranges input file for inventor name clustering
	%% in the next stage.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	clear;
	fid=fopen('./_attribute_matrix','r');
	temp=textscan(fid,'%d %d %d','delimiter',',','bufsize',2999999);
	fclose(fid);
	X=sparse(double(temp{1}),double(temp{2}),double(temp{3}));
	clear temp;
	fid=fopen('./_attribute_metric','r');
	temp=textscan(fid,'%s %f %s %s %s %s %s','delimiter','\t','bufsize',2999999);
	fclose(fid);
	patno=temp{1};
	inventor=temp{3};
	assignee=temp{4};
	class=temp{5};
	lastname=temp{6};
	inventor_id=temp{7};
	clear temp;

	XX=X>0;
	XXt=XX';

	fid=fopen(sprintf('./_attribute_dictionary'),'r');
	words=textscan(fid,'%*d %s %*d','delimiter','\t','bufsize',9999999);
	fclose(fid);
	words=words{1};

	fid=fopen(sprintf('./_disambiguator_input.csv'),'r');
	temp=textscan(fid,'%s %s %s %s %*s %*s %*s %*s %*s %*s %*s %*s','delimiter','\t','bufsize',9999999);
	fclose(fid);
	key=temp{1};
	name1=temp{2};
	name2=temp{3};
	name3=temp{4};
	namex=name1;
	for i=1:length(name1)
		namex{i}=[name1{i} ' ' name2{i} ' ' name3{i}];
	end
	namex=regexprep(namex,'[\ ]{1,}',' ');
	namex=strtrim(namex);
	clear temp i;



