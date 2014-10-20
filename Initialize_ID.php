<?php
	ini_set('memory_limit','20000M');
	$input="_disambiguator_input.csv";

	$outf1=fopen("_attribute_line","w");
	$outf2=fopen("_attribute_metric","w");

	$inf=fopen($input,"r");
	$patno=Array();
	for ($n=0;($line=fgets($inf))!=false;$n++)
	{
		$linex=explode("\t",$line);
		$patno[]=$linex[4];
	}
	fclose($inf);

	$inf=fopen($input,"r");
	$inventors=[];
	$inventor_id=[];
	$inventors_start=[];
	$inv="";
	$start=0;
	for($n=0;$n<sizeof($patno)-1;)
	{
		for(;$n<sizeof($patno)+1;$n++)
		{
			if($n==sizeof($patno)||strcmp($patno[$start],$patno[$n]))
			{
				$end=$n-1;
				$inv=trim($inv);
				$inventors[$start]=$inv;;
				$inv="";
				$start=$n;
				if($n<sizeof($patno))
					$n--;
				else
					break;
			}
			else
			{
				$inventors_start[]=$start;
				$inventor_id[]=$patno[$n]."-".($n-$start+1);
				$line=fgets($inf);
				$linex=explode("\t",$line);
				$invx=$linex[1]." ".$linex[2]." ".$linex[3];
				$invx=preg_replace("/[\t\r\n]/"," ",$invx);
				$invx=preg_replace("/[\ ]+/","_",$invx);
				$invx=trim($invx);
				$inv=$inv.$invx." ";
			}
		}
	}

	fclose($inf);
	$inf=fopen($input,"r");

	for ($n=0;($line=fgets($inf))!=false;$n++)
	{

		$line=str_replace("\n","",$line);
		$line=explode("\t",$line);
		$name=$line[1]." ".$line[2]." ".$line[3];
		$name=preg_replace("/[\t\r\n]/"," ",$name);
		$name=preg_replace("/[\ ]+/","_",$name);
		$assignee=$line[11];
		$assignee=preg_replace("/[\t\r\n]/"," ",$assignee);

		$out2=sprintf("%s\t0\t%s\t%s\t%s\t%s\t%s",$line[4],$name,$assignee,$line[5],$line[3],$inventor_id[$n]);
		$out2=trim($out2);
		$out2=strtoupper($out2);
		fprintf($outf2,"%s\n",$out2);

		$class=$line[5];
		$class=strtolower($class);

		$assignee=$line[11];
		$city=$line[7];

		$lastname=preg_replace("/[\t\r\n]/"," ",$line[3]);
		$lastname=preg_replace("/[\ ]+/","_",$lastname);
		$lastname=strtolower($lastname);
		$assignee=preg_replace("/[\t\r\n]/"," ",$line[11]);
		$assignee=preg_replace("/[\ ]+/","_",$assignee);
		$assignee=strtolower($assignee);
		$city=preg_replace("/[\t\r\n]/"," ",$line[7]);
		$city=preg_replace("/[\ ]+/","_",$city);
		$city=strtolower($city);

		$name=strtolower($name);

		if(strlen($name)>0)
			$firstname_letter=substr($name,0,1);
		else
			$firstname_letter="";

		if(strlen($assignee)>0 && strcmp($assignee,"|"))
			$out1=sprintf("FN%s LN%s NM%s AS%s CT%s CL%s",$firstname_letter,$lastname,$name,$assignee,$city,$class);
		else
			$out1=sprintf("FN%s LN%s NM%s AS CT%s CL%s",$firstname_letter,$lastname,$name,$city,$class);

		$out1=preg_replace("/[\ ]+/"," ",$out1);
		$out1=trim($out1);
		fprintf($outf1,"%s\n",$out1);
	}
	fclose($inf);
	fclose($outf1);
	fclose($outf2);

?>
