
start_time=now;

%%% matrix and sub-matrix transpose to expedite memory access %%%
%%% matrix and sub-matrix transpose to expedite memory access %%%

keep=strncmp(words,'LN',2);
XXLN=XX(:,keep);
XXLNt=XXLN';

keep=strncmp(words,'FN',2);
XXFN=XX(:,keep);
XXFNt=XXFN';

keep=strncmp(words,'NM',2);
XXNM=XX(:,keep);
XXNMt=XXNM';

keep=strncmp(words,'AS',2);
XXAS=XX(:,keep);
XXASt=XXAS';

keep=strncmp(words,'CT',2);
XXCT=XX(:,keep);
XXCTt=XXCT';

keep=strncmp(words,'CL',2);
XXCL=XX(:,keep);
XXCLt=XXCL';

patno_1=1:length(patno);
patno_1=patno_1';

%%% matrix and sub-matrix transpose to expedite memory access %%%
%%% matrix and sub-matrix transpose to expedite memory access %%%

step=ones(length(patno),1)>0;

inventor_id_new=inventor_id;
namey=namex;

for index=1:length(patno)
	if ~step(index)
		continue;
	end

		target=XXLNt(:,index);
		C_lastname=XXLN(:,target);

		target=XXFNt(:,index);
		C_firstname=XXFN(:,target);

		target=XXNMt(:,index);
		C_name=XXNM(:,target);

		target=XXASt(:,index);
		C_assignee=XXAS(:,target);

		target=XXCTt(:,index);
		C_city=XXCT(:,target);

		target=XXCLt(:,index);
		C_class=XXCL(:,target);

		lump_index_2=step&C_name&((C_assignee|C_city|C_class)); %% perfect match %%
		lump_index_1=step&C_firstname&C_lastname&((C_assignee&C_city&C_class)); %% imperfect match %%
		lump_patno_2=patno(lump_index_2);
		lump_patno_1=patno(lump_index_1);
		lump_index_2=find(lump_index_2);
		lump_index_1=find(lump_index_1);
		lump_index_1_=[];

%% if exists in patno_1 and in patno_2, then take out patno_1 (full match)
		for ii=1:length(lump_index_1)
			if ~any(strcmp(lump_patno_1{ii},lump_patno_2))
				lump_index_1_(end+1)=lump_index_1(ii);
			end
		end
%% if exists in patno_1 and in patno_2, then take out patno_1 (full match)

		lump_index=union(lump_index_1_,lump_index_2);

%%% find repeated patent NO.s adopted %%%

		lump_initial=lump_index;

		lump_after=lump_initial;

	if length(lump_initial)>0

		for indexx=1:length(lump_initial)

		indexy=lump_initial(indexx);

		target=XXLNt(:,indexy);
		C_lastname=XXLN(:,target);

		target=XXFNt(:,indexy);
		C_firstname=XXFN(:,target);

		target=XXNMt(:,indexy);
		C_name=XXNM(:,target);

		target=XXASt(:,indexy);
		C_assignee=XXAS(:,target);

		target=XXCTt(:,indexy);
		C_city=XXCT(:,target);

		target=XXCLt(:,indexy);
		C_class=XXCL(:,target);

		lump_index_2=step&C_name&((C_assignee|C_city|C_class)); %% perfect match %%
		lump_index_1=step&C_firstname&C_lastname&((C_assignee&C_city&C_class)); %% imperfect match %%
		lump_patno_2=patno(lump_index_2);
		lump_patno_1=patno(lump_index_1);
		lump_index_2=find(lump_index_2);
		lump_index_1=find(lump_index_1);
		lump_index_1_=[];

	%% if exists in patno_1 and in patno_2, then take out patno_1 (full match)
		for ii=1:length(lump_index_1)
			if ~any(strcmp(lump_patno_1{ii},lump_patno_2))
				lump_index_1_(end+1)=lump_index_1(ii);
			end
		end
	%% if exists in patno_1 and in patno_2, then take out patno_1 (full match)

		lump_index=union(lump_index_1_,lump_index_2);

		lump_after_2=lump_index;

		lump_after=union(lump_after,lump_after_2);

		end

	end

	lump=union(lump_initial,lump_after);

	if length(lump)==0
		step(index)=false;
		continue;
	end

	for i=length(lump):-1:1
		step(lump(i))=false;
		inventor_id_new{lump(i)}=inventor_id_new{index};
		namey{lump(i)}=namey{index};
	end
end

fprintf('\n\n');

	fid0=fopen('_disambiguator_output.tsv','w');
	for i=1:length(inventor_id_new)
		fprintf(fid0,'%s\t%s\n',key{i},inventor_id_new{i});
	end
	fclose(fid0);



