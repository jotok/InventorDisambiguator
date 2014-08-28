<?php

// Removed all code that was trying to use multiple threads since the PCTNL library
// is not available for Windows

function scan_dir_file($path)
{
	$dirr=scandir($path);
	unset($dirr[0]);
	unset($dirr[1]);
	$dirr=array_merge($dirr);
	for($i=0;$i<sizeof($dirr);$i++)
		$dirr[$i]=$path."\\".$dirr[$i];
	return $dirr;
}

function combine_newsdict($file1,$words1,$words2)
{
	if(strlen($words1)>0)
	{
		$words1=substr($words1,0,strlen($words1)-1);
		$words1=explode("\n",$words1);
	}
	else
		$words1=Array();

	if(strlen($words2)>0)
	{
		$words2=substr($words2,0,strlen($words2)-1);
		$words2=explode("\n",$words2);
	}
	else
		$words2=Array();

	$dict=array_unique(array_merge($words1,$words2));
	$dict=array_merge($dict);

	$outf=fopen($file1,"w");
	for ($i=0;$i<sizeof($dict);$i++)
		fprintf($outf,"%s\n",$dict[$i]);
	fclose($outf);

	return;
}

function combine_threaded_newsdict_multi_thread($path,$outfn)
{
	$dirr=scan_dir_file($path);
	$pid_arr=Array();

	$layer=-1;
	$total=0;

	while(sizeof($dirr)+sizeof($pid_arr)>1)
	{
		$layer++;

		while(sizeof($dirr)>1)
		{
		
			$file1=$dirr[0];
			$file2=$dirr[1];
			$words1=file_get_contents($dirr[0]);
			$words2=file_get_contents($dirr[1]);
			shell_exec("del /q ".$dirr[0]);
			shell_exec("del /q ".$dirr[1]);
			unset($dirr[0]);
			unset($dirr[1]);
			$dirr=array_merge($dirr);

/*			$pid=pcntl_fork();
			if ($pid==-1)
				die('could not fork');
			else
			{
				if ($pid)
				{
					$pid_arr[]=$pid;
				}
				else
				{
					combine_newsdict($file1,$words1,$words2);
					exit(0);
				}
			} */
			combine_newsdict($file1,$words1,$words2);
		}
/*		while(count($pid_arr) > 0)
		{
			$myId = pcntl_waitpid(-1, $status, WNOHANG);
			foreach($pid_arr as $key => $pid)
			{
				if($myId==$pid)
				{
					unset($pid_arr[$key]);
					$dirr=scan_dir_file($path);
					break;
				}
			}
			usleep(100);
		} */
	}

/*	while(count($pid_arr) > 0)
	{
		$myId = pcntl_waitpid(-1, $status, WNOHANG);
		foreach($pid_arr as $key => $pid)
			if($myId==$pid)
				unset($pid_arr[$key]);
		usleep(100);
	} */

	$dict=file_get_contents($dirr[0]);
	$dict=substr($dict,0,strlen($dict)-1);
	$dict=explode("\n",$dict);

	$dictt=Array();

	$outf=fopen($outfn,"w");
	for($i=0;$i<sizeof($dict);$i++)
	{
		fprintf($outf,"%d\t%s\t%d\n",$i+1,$dict[$i],"1");
		$dictt[$dict[$i]]=$i+1;
	}
	fclose($outf);

	return $dictt;

}

function combine_threaded_newspara($path,$b,$outfn)
{
	$outf=fopen($outfn,"w");
	for ($i=0;$i<$b;$i++)
	{
		$inf=fopen($path."\\".$i,"r");
		for (;($content=fgets($inf))!=false;)
			fprintf($outf,"%s",$content);
		fclose($inf);
	}
	fclose($outf);

	return;
}

function generate_newspara_multi_thread($path,$bi,$infn,$start,$end,$dict)
{

	$inf=fopen($infn,"r");
	$outf=fopen($path."\\".$bi,"w");

	for ($n=0;($content=fgets($inf))!=false;$n++)
	{
		if ($n<$start or $n>$end)
			continue;

		$content=str_replace("\n","",$content);
		$content=explode(" ",$content);
		$content=array_count_values($content);
		$words=array_keys($content);

		foreach($words as $word)
			fprintf($outf,"%d,%s,%d\n",$n+1,$dict[$word],$content[$word]);
	}
	fclose($inf);
	fclose($outf);


	return;
}

function generate_newsdict_multi_thread($path,$bi,$infn,$start,$end)
{
	$inf=fopen($infn,"r");
	$outf=fopen($path."\\".$bi,"w");

	$method="load_2";

	$dict=Array();

	if(!strcmp($method,"load_1"))
	{
		for ($n=0;($content=fgets($inf))!=false;$n++)
		{
			if ($n<$start or $n>$end)
				continue;
			$content=str_replace("\n","",$content);
			$content=explode(" ",$content);
			$dict=array_unique(array_merge($dict,$content));
		}
	}
	if(!strcmp($method,"load_2"))
	{
		$temp="";
		for ($n=0;($content=fgets($inf))!=false;$n++)
		{
			if ($n<$start or $n>$end)
				continue;
			$content=str_replace("\n","",$content);
			$temp.=$content." ";
		}
		$temp=substr($temp,0,strlen($temp)-1);
		$temp=explode(" ",$temp);
		$dict=array_unique(array_merge($temp));
	}

	fclose($inf);

	$dict=array_merge($dict);
	$dictt=Array();

	for ($i=0;$i<sizeof($dict);$i++)
		fprintf($outf,"%s\n",$dict[$i]);
	fclose($outf);

	return;
}

function generate_newsdict($infn,$outfn,$b)
{
	$inf=fopen($infn,"r");
	for ($n=0;($line=fgets($inf))!=false;)
		$n++;
	fclose($inf);

	$intervals=Array();
	$block_size_base=floor($n/$b);
	$block_residual=$n%$b;
	for ($i=0;$i<$b;$i++)
	{
		$block_size=$block_size_base+($i<$block_residual);
		if ($i<=$block_residual)
			$block_start=$i*ceil($n/$b);
		else
			$block_start=($i-$block_residual)*$block_size_base+$block_residual*ceil($n/$b);
		$block_end=$block_start+$block_size-1;
		$intervals[$i]["start"]=$block_start;
		$intervals[$i]["end"]=$block_end;
	}

	$path=".\\tmp\\php".getmypid();
	shell_exec("md ".$path);

	$starttime=microtime(TRUE);
	$elapsed = microtime(TRUE) - $starttime;

	$pid_arr=Array();
	for ($i=0;$i<$b;$i++)
	{
/*		$pid=pcntl_fork();
		if ($pid==-1)
			die('could not fork');
		else
		{
			if ($pid)
				$pid_arr[$i]=$pid;
			else
			{
				generate_newsdict_multi_thread($path,$i,$infn,$intervals[$i]["start"],$intervals[$i]["end"]);
				exit(0);
			}
		} */
		generate_newsdict_multi_thread($path,$i,$infn,$intervals[$i]["start"],$intervals[$i]["end"]);
	}

/*	while(count($pid_arr) > 0)
	{
		$myId = pcntl_waitpid(-1, $status, WNOHANG);
		foreach($pid_arr as $key => $pid)
		{
			if($myId==$pid)
				unset($pid_arr[$key]);
		}
		usleep(100);
	} */

	$dict=combine_threaded_newsdict_multi_thread($path,$outfn);

	shell_exec("rmdir /s /q ".$path);

	return $dict;
}

function generate_newspara($infn,$outfn,$dict,$b)
{
	$inf=fopen($infn,"r");
	for ($n=0;($line=fgets($inf))!=false;)
		$n++;
	fclose($inf);

	$intervals=Array();
	$block_size_base=floor($n/$b);
	$block_residual=$n%$b;
	for ($i=0;$i<$b;$i++)
	{
		$block_size=$block_size_base+($i<$block_residual);
		if ($i<=$block_residual)
			$block_start=$i*ceil($n/$b);
		else
			$block_start=($i-$block_residual)*$block_size_base+$block_residual*ceil($n/$b);
		$block_end=$block_start+$block_size-1;
		$intervals[$i]["start"]=$block_start;
		$intervals[$i]["end"]=$block_end;
	}

	$path=".\\tmp\\php".getmypid();
	shell_exec("md ".$path);

	$starttime=microtime(TRUE);
	$elapsed = microtime(TRUE) - $starttime;

	$pid_arr=Array();
	for ($i=0;$i<$b;$i++)
	{
/*		$pid=pcntl_fork();
		if ($pid==-1)
			die('could not fork');
		else
		{
			if ($pid)
				$pid_arr[$i]=$pid;
			else
			{
				generate_newspara_multi_thread($path,$i,$infn,$intervals[$i]["start"],$intervals[$i]["end"],$dict);
				exit(0);
			}
		} */
		generate_newspara_multi_thread($path,$i,$infn,$intervals[$i]["start"],$intervals[$i]["end"],$dict);
	}

	while(count($pid_arr) > 0)
	{
		$myId = pcntl_waitpid(-1, $status, WNOHANG);
		foreach($pid_arr as $key => $pid)
		{
			if($myId==$pid)
				unset($pid_arr[$key]);
		}
		usleep(100);
	}

	combine_threaded_newspara($path,$b,$outfn);

	shell_exec("rmdir /s /q ".$path);

	return 0;
}

	ini_set('memory_limit','12000M');
	date_default_timezone_set('America/Los_Angeles');

	$starttime=microtime(TRUE);

	$number_of_threads_newsdict=1;
	$number_of_threads_newspara=1;

	$dict=generate_newsdict("_attribute_line","_attribute_dictionary",$number_of_threads_newsdict);

	generate_newspara("_attribute_line","_attribute_matrix",$dict,$number_of_threads_newspara);

	$elapsed = microtime(TRUE) - $starttime;

?>
