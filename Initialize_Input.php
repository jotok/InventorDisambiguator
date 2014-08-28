<?php

	$input=$argv[1];
	$outf=fopen("_disambiguator_input.csv","w");

	$inf=fopen($input,"r");
	$patno=Array();
	for ($n=0;($line=fgets($inf))!=false;$n++)
	{
		$linex=explode("\t",$line);
		if(strcmp($linex[2],"1"))
			fprintf($outf,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
				$linex[0],
				$linex[3],
				$linex[4],
				$linex[5],
				$linex[6],
				$linex[7],
				$linex[8],
				$linex[9],
				$linex[10],
				$linex[11],
				$linex[12],
				$linex[13]
			);
	}
	fclose($inf);
	fclose($outf);

?>
