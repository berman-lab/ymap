<?php
	session_start();
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<HTML>
<HEAD>
	<style type="text/css">
		body {font-family: arial;}
	</style>
	<script>	
		// will return the number of checked draw check boxes
		function getCheckedNumber(checkBoxes) {
			// getting the checkboxes
			var i = 0; // loop variable
			var count = 0; // will count the number of checked checkboxes
			for (i = 0; i < checkBoxes.length; i++) {
				if (checkBoxes[i].checked) {
					count += 1;				
				}
			}
			return count;
		}

		function chromosomCheck() { 
			// getting the checkboxes
			var checkBoxes = document.getElementsByClassName("draw");
			var count = getCheckedNumber(checkBoxes);
			var i = 0; // loop variable		
			// disable all unchecked boxes if we have reached limit			
			if (count == 50) {
				for (i = 0; i < checkBoxes.length; i++) {
					if (!checkBoxes[i].checked) {
						checkBoxes[i].disabled = true;				
					}
				}
			}
			else { // enable checkboxes
				for (var i = 0; i < checkBoxes.length; i++) {
					checkBoxes[i].disabled = false;				
				}	
			}
		}
	</script>
<meta http-equiv="content-type" content="text/html; charset=iso-8859-1">
<title>Install genome into pipeline.</title>
</HEAD>
<?php
// current directory
	require_once '../constants.php';
	include_once 'process_input_files.genome.php';

	$fileName = filter_input(INPUT_POST, "fileName", FILTER_SANITIZE_STRING);
	$user     = filter_input(INPUT_POST, "user",     FILTER_SANITIZE_STRING);
	$genome   = filter_input(INPUT_POST, "genome",   FILTER_SANITIZE_STRING);
	$key      = filter_input(INPUT_POST, "key",      FILTER_SANITIZE_STRING);

//// Debugging
// echo getcwd() . "\n";
// echo $fileName."\n";
// echo $user."\n";
// echo $genome."\n";
// echo $key."\n";

	$fasta_name = $genome.".fasta";
	if (($fileName == ".") || ($fileName == "..") || ($user == ".") || ($user == "..") || ($genome == ".") || ($genome == "..") || ($key == ".") || ($key == "..")) {
		echo "Invalid input data";
		exit;
	}

// Initialize 'process_log.txt' file.
	$logOutputName = "../users/".$user."/genomes/".$genome."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'w');
	fwrite($logOutput, "Log file initialized\n");
	fwrite($logOutput, "Running 'scripts_genomes/genome.install_1.php'.\n");

// Initialize 'condensed_log.txt' file.
	$condensedLogOutputName = "../users/".$user."/genomes/".$genome."/condensed_log.txt";
	$condensedLogOutput     = fopen($condensedLogOutputName, 'w');
	fwrite($condensedLogOutput, "Initializing.\n");
	fclose($condensedLogOutput);
	chmod($condensedLogOutputName,0744);


	// Generate 'reference.txt' file containing:
	//      one line; file name of reference FASTA file.
	fwrite($logOutput, "\tGenerating 'reference.txt' file.\n");
	$outputName       = "../users/".$user."/genomes/".$genome."/reference.txt";
	if (file_exists($outputName)) {
		$fileContents = file_get_contents($outputName);
		unlink($outputName);
		$output       = fopen($outputName, 'w');
		fwrite($output, $fasta_name);
	} else {
		$output       = fopen($outputName, 'w');
		fwrite($output, $fasta_name);
	}
	fclose($output);
	unset($outputName);
	unset($output);

	// Generate 'name.txt' file containing:
	//		one line; name of genome.
	fwrite($logOutput, "\tGenerating 'name.txt' file.\n");
	$outputName       = "../users/".$user."/genomes/".$genome."/name.txt";
	if (file_exists($outputName)) {
		$fileContents = file_get_contents($outputName);
		unlink($outputName);
		$output       = fopen($outputName, 'w');
		fwrite($output, $fileContents);
	} else {
		$output       = fopen($outputName, 'w');
		fwrite($output, str_replace("_"," ",$genome));
	}
	fclose($output);
	unset($outputName);
	unset($output);

	// Generate 'upload_size.txt' file to contain the size of the uploaded file (irrespective of format) for display in "Manage Datasets" tab.
	$genomePath      = "../users/".$user."/genomes/".$genome."/";
	$outputName      = $genomePath."upload_size_1.txt";
	$output          = fopen($outputName, 'w');
	$fileSizeString  = filesize($genomePath.$fileName);
	fwrite($output, $fileSizeString);
	fclose($output);
	chmod($outputName,0755);
	fwrite($logOutput, "\tGenerated 'upload_size_1.txt' file.\n");

		// Process uploaded file.
		$name        = str_replace("\\", ",", $fileName);
		rename($genomePath.$name,$genomePath.strtolower($name));
		$name        = strtolower($name);
		$ext         = strtolower(pathinfo($name, PATHINFO_EXTENSION));
		$filename    = strtolower(pathinfo($name, PATHINFO_FILENAME));
		fwrite($logOutput, "\tKey         : ".$key."\n");
		fwrite($logOutput, "\tDatafile    : '".$name."'.\n");
		fwrite($logOutput, "\tFilename    : '".$filename."'.\n");
		fwrite($logOutput, "\tExtension   : '".$ext."'.\n");
		fwrite($logOutput, "\tScript_path : '".getcwd()."'.\n");
		fwrite($logOutput, "\tPath        : '".$genomePath."'.\n");

		// Generate 'upload_size.txt' file to contain the size of the uploaded file (irrespective of format) for display in "Manage Datasets" tab.
		$output2Name    = $genomePath."upload_size_1.txt";
		$output2        = fopen($output2Name, 'w');
		$fileSizeString = filesize($genomePath.$name);
		fwrite($output2, $fileSizeString);
		fclose($output2);
		chmod($output2Name,0755);
		fwrite($logOutput, "\tGenerated 'upload_size_1.txt' file.\n");

		// Process the uploaded file.
		process_input_files_genome($ext,$name,$genomePath,$key,$user,$genome,$output, $condensedLogOutput,$logOutput, $fasta_name);
		$fileName = $fasta_name;

	$file_path  = "../users/".$user."/genomes/".$genome."/".$fileName;
	// Process FASTA file for chromosome count, names, and lengths.
	fwrite($logOutput, "\tReading chromosome count, names, and lengths from FASTA.\n");
	$chr_count   = 0;
	$chr_names   = array();
	$chr_lengths = array();	
	$chr_length = 0;
	$fileHandle = fopen($file_path, "r") or die("Couldn't open fasta file after first reformat");
	$firstLine = true;	
	if ($fileHandle) {
		while (!feof($fileHandle)) {
			$buffer = fgets($fileHandle, 4096);
			if ($buffer[0] == '>') {
				if ($chr_length != 0) {
					fwrite($logOutput, "\t" . 'finished processing chromosome ' . $chr_name . ' length: ' . $chr_length . PHP_EOL);
				}
				// chromosome name is the header string (starting with ">"), after trimming trailing whitespace characters.				
				$line_parts = explode(" ",$buffer);
				$chr_name   = str_replace(array("\r","\n"),"",$line_parts[0]);
				$chr_name   = substr($chr_name,1,strlen($chr_name)-1);
				array_push($chr_names, $chr_name);
				$chr_count  += 1;
				// pushing chromosome length only if it's not the first line to avoid pushing 0
				if ($firstLine == false) {
					// pushing latest chromosome length and reseting count
					array_push($chr_lengths, $chr_length);
					$chr_length = 0;
				} 
				else {
					$firstLine = false;
				}		
			}
			else {
				// adding char count without line end or whitespaces at start
				$chr_length += strlen(trim($buffer));	
			}
		}
		// adding last chr_length, pushing zero if it doesn't exist
		if ($chr_length > 0)
		{
			array_push($chr_lengths, $chr_length);
			fwrite($logOutput, "\t" . 'finished processing chromosome ' . $chr_name . ' length: ' . $chr_length . PHP_EOL);
		}
		fclose($handle);
	}
	unset($line);
	unset($chr_length);
	unset($line_parts);
	unset($chr_name);
	fwrite($logOutput, "\tnum chrs    : '$chr_count'.\n");

	// Store variables of interest into $_SESSION.
	fwrite($logOutput, "\tStoring PHP session variables.\n");
	$_SESSION['genome_'.$key]      = $genome;
	$_SESSION['fileName_'.$key]    = $fileName;
	$_SESSION['chr_count_'.$key]   = $chr_count;
	// saving all chr details in files to avoid overloading $_SESSION
	file_put_contents("../users/".$user."/genomes/".$genome."/chr_names.json",json_encode($chr_names));
	file_put_contents("../users/".$user."/genomes/".$genome."/chr_lengths.json",json_encode($chr_lengths));

//// Debugging output of all variables.
//	print_r($_SESSION);
//	print_r($GLOBALS);

	// The following section defines a form for collecting the information needed to build the last of the genome setup files.
	fwrite($logOutput, "\tGenerating form to request centromere location and other genome specific data from the user.\n");
?>

<BODY onload = "parent.parent.resize_genome('<?php echo $key; ?>', 150)<?php
//	$genomePath      = "../users/".$user."/genomes/".$genome."/";
//	$outputName      = $genomePath."upload_size_1.txt";
//	$output          = fopen($outputName, 'w');
//	$fileSizeString  = filesize($genomePath.$fileName);

	$sizeFile_1   = "upload_size_1.txt";
	$handle       = fopen($sizeFile_1,'r');
	$sizeString_1 = trim(fgets($handle));
	fclose($handle);
	if ($sizeString_1 !== "") {
		echo "; parent.parent.update_genome_file_size('".$key."','".$sizeString_1."');";
	}
?>">
<?php   
	echo "<font color=\"red\" size=\"2\">Fill in genome details:</font><br>";
	echo "<font color=\"black\" size=\"2\">Label can be up to 6 characters</font><br>";
	if ($chr_count > 50 && $chr_count <= 300) { 
		echo "<font color=\"black\" size=\"2\">Ymap can only display up to 50 scaffolds.<br>The longest 50 have been automatically selected, though you can change the selection.</font>";
	}
	else if ($chr_count > 300) {
		echo "<font color=\"black\" size=\"2\">Ymap can only work with up to 300 scaffolds, and visualize only up to 50. The reference you uploaded had " . $chr_count ." scaffolds. 						The longest 300 are available for choosing and the longest 50 have been automatically selected, though you can change the selection in the form below.</font>";
	}
	echo "<form name=\"chromSelect\" action=\"genome.install_2.php\" method=\"post\">";	
		echo "<table border=\"0\">";
		echo "<tr>";
			echo "<th rowspan=\"2\"><font size=\"2\">Use</font></th>";
			echo "<th rowspan=\"2\"><font size=\"2\">FASTA entry name</font></th>";
			echo "<th rowspan=\"2\"><font size=\"2\">Label</font></th>";
			echo "<th colspan=\"2\"><font size=\"2\">Centromere</font></th>";
			echo "<th rowspan=\"2\"><font size=\"2\">rDNA</font></th>";
			echo "<th rowspan=\"2\"><font size=\"2\">Size(BP)</font></th>";
		echo "</tr>";
		echo "<tr>";
			echo "<th><font size=\"2\">start bp</font></th>";
			echo "<th><font size=\"2\">end bp</font></th>";
		echo "</tr>";

	// if the chr_count is above 50 sorting and getting the size of the 50 chromosome to use as a reference whether to check or uncheck chromosomes
	// also getting the size of the 			
	if ($chr_count > 50) {
		$chr_lengthsTemp = $chr_lengths;		
		rsort($chr_lengthsTemp);
		// getting cutoff value for the 50 longest chromsomes
		$lowestSize = $chr_lengthsTemp[49];
		if ($chr_count > 300) {				
			$lowestSizeDisplay = $chr_lengthsTemp[299];
		}
		// clearing the copied array
		unset($chr_lengthsTemp);
		$countUsed = 50; // will decrease for each checked chromosome
	}
	for ($chr=0; $chr<$chr_count; $chr+=1) {
		$chrID = $chr+1;
		// in case we have to many chromsomes dispaying only the 900 longest ones, so jumping lower size chromosomes			
		if ($chr_count > 300 && $chr_lengths[$chr] < $lowestSizeDisplay) {
			continue;				
		}
		echo "\t\t<tr>\n";
		// disabling chromosomes with length of 0 from been checked
		if ($chr_lengths[$chr] == 0) {
			echo "\t\t\t<td align=\"middle\"><input type=\"checkbox\" class=\"draw_disabled\" name=\"draw_{$chrID}\" disabled></td>\n";
		}
		// checking all if number of chromsomes is less or equal to 50
		else if ($chr_count <= 50) {
			echo "\t\t\t<td align=\"middle\"><input type=\"checkbox\" class=\"draw\" name=\"draw_{$chrID}\" checked></td>\n";
		}
		else {
			// checking only the 50 longest
			if ($chr_lengths[$chr] >= $lowestSize && $countUsed >= 1) {
				echo "\t\t\t<td align=\"middle\"><input type=\"checkbox\" class=\"draw\" name=\"draw_{$chrID}\" onchange=\"chromosomCheck()\" checked></td>\n";
				// reduce the count of used to avoid checking more then 50 if multiple have the same size
				$countUsed -= 1;
			}
			else {
				echo "\t\t\t<td align=\"middle\"><input type=\"checkbox\" class=\"draw\" onchange=\"chromosomCheck()\" name=\"draw_{$chrID}\" disabled></td>\n";
			}
		}
		echo "\t\t\t<td><font size=\"2\">{$chr_names[$chr]}</font></td>\n";
		echo "\t\t\t<td><input type=\"text\"     name=\"short_{$chrID}\"    value=\"Chr{$chrID}\" size=\"6\"maxlength=\"6\" ></td>\n";
		echo "\t\t\t<td><input type=\"text\"     name=\"cenStart_{$chrID}\" value=\"0\"           size=\"6\"></td>\n";
		echo "\t\t\t<td><input type=\"text\"     name=\"cenEnd_{$chrID}\"   value=\"0\"           size=\"6\"></td>\n";
		echo "\t\t\t<td align=\"middle\" ><input type=\"radio\"    name=\"rDNAchr\"      value=\"{$chrID}\"></td>\n";
		echo "\t\t\t<td><font size=\"2\">{$chr_lengths[$chr]}</font></td>\n";
		echo "\t\t</tr>\n";
	}

	echo "</table><br>";
	echo "<font size=\"2\">";
	echo "Ploidy = <input type=\"text\" name=\"ploidy\" value=\"2.0\" size=\"6\"><br>";
	echo "rDNA (start = <input type=\"text\" name=\"rDNAstart\" value=\"0\" size=\"6\">; end = <input type=\"text\" name=\"rDNAend\" value=\"0\" size=\"6\">)<br>";
	echo "Futher annotations to add to the genome? <input type=\"text\" name=\"annotation_count\" value=\"0\" size=\"6\"><br>";
	echo "Select <input type=\"checkbox\" name=\"expression_regions\"> if a tab-delimited-text file listing ORF coordinates is available.<br>";
	echo "</font>";
	echo "<br>";
	echo "<input type=\"submit\" value=\"Save genome details...\">";
	echo "<input type=\"hidden\" id=\"key\" name=\"key\" value=\"".  $key . "\">";
	echo "</form>"; 
?>

</BODY>
</HTML>
<?php
	fwrite($logOutput, "\t'scripts_genomes/genome.install_1.php' has completed.\n");
	fclose($logOutput);
?>
