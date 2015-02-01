<?php
	session_start();
	if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script> <?php } else { $user = $_SESSION['user']; }
	require_once 'constants.php';
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";
	error_reporting(E_ALL);
	ini_set('display_errors', 1);

	$bad_chars       = array(".", ",", "\\", "/", " ");
	$hapmap          = str_replace($bad_chars,"_",trim( filter_input(INPUT_POST, "hapmap", FILTER_SANITIZE_STRING) ));
	$user            = $_SESSION['user'];
	$genome          =        filter_input(INPUT_POST, "genome",             FILTER_SANITIZE_STRING);
	$referencePloidy = (float)filter_input(INPUT_POST, "referencePloidy",    FILTER_SANITIZE_NUMBER_FLOAT);
	$project1        =        filter_input(INPUT_POST, "project1",           FILTER_SANITIZE_STRING);
	$project2        =        filter_input(INPUT_POST, "project2",           FILTER_SANITIZE_STRING);
	$colorA          =        filter_input(INPUT_POST, "homolog_a_color",    FILTER_SANITIZE_STRING);
	$colorB          =        filter_input(INPUT_POST, "homolog_b_color",    FILTER_SANITIZE_STRING);
	if ($referencePloidy == 2) {
		$hapmap_description = filter_input(INPUT_POST, "hapmap_description", FILTER_SANITIZE_STRING);
	}

	$dir1      = "../users/".$user."/hapmaps";
	$dir2      = "../users/".$user."/hapmaps/".$hapmap;
	$dir3      = "../users/default/hapmaps/".$hapmap;

	$currentPath = getcwd();

// Create the hapmap folder inside the user's hapmaps directory
	if (file_exists($dir2)) {
	} else {
		mkdir($dir2);
		chmod($dir2,0777);
	}

// Initialize 'process_log.txt' file.
	$logOutputName = "../users/".$user."/hapmaps/".$hapmap."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Running 'php/hapmap.install_3.php'.\n");
	fwrite($logOutput, "\tuser            = ".$user."\n");
	fwrite($logOutput, "\thapmap          = ".$hapmap."\n");
	fwrite($logOutput, "\tgenome          = ".$genome."\n");
	fwrite($logOutput, "\treferencePloidy = ".$referencePloidy."\n");
	fwrite($logOutput, "\tproject1        = ".$project1."\n");
	fwrite($logOutput, "\tproject2        = ".$project2."\n");
	fwrite($logOutput, "\tcolorA          = ".$colorA."\n");
	fwrite($logOutput, "\tcolorB          = ".$colorB."\n");

// Create 'colors.txt' file to contain colors used in haplotype figures.
	$handleName = "../users/".$user."/hapmaps/".$hapmap."/colors.txt";
	if (file_exists($handleName)) {
	} else {
		$handle     = fopen($handleName, 'w');
		fwrite($handle, $colorA."\n".$colorB);
		fclose($handle);
	}

// Create 'genome.txt' file to contain genome used in haplotype.
	$handleName = "../users/".$user."/hapmaps/".$hapmap."/genome.txt";
	if (file_exists($handleName)) {
	} else {
		$handle     = fopen($handleName, 'w');
		fwrite($handle, $genome);
		fclose($handle);
	}

// Create 'parent.txt' file to contain parent genome used in haplotype.
	$handleName = "../users/".$user."/hapmaps/".$hapmap."/parent.txt";
	if (file_exists($handleName)) {
	} else {
		$handle     = fopen($handleName, 'w');
		if ($referencePloidy == 2) {
			fwrite($handle, $project1);
		} else {
			fwrite($handle, $project1."\n".$project2);
		}
		fclose($handle);
	}

// Initialize 'haplotype.txt' file to hold haplotype entry descriptions.
	if ($referencePloidy == 2) {
		$haplotypeFileName = "../users/".$user."/hapmaps/".$hapmap."/haplotypeMap.txt";
		$haplotypeFile     = fopen($haplotypeFileName, 'a');
		fwrite($haplotypeFile, $project1."\n".$project2."\n".$hapmap_description."\n");
		fclose($haplotypeFile);
	} else {
		// for haplotype map derived from two haploid strains, there will not be more entries after construction.
		// Because of this, 'parent.txt' will contain sufficient information to define haplotype map.
	}

// Generate 'working.txt' file to let pipeline know that processing is underway.
	$handleName      = "../users/".$user."/hapmaps/".$hapmap."/working.txt";
	$handle          = fopen($handleName, 'w');
	$startTimeString = date("Y-m-d H:i:s");
	fwrite($handle, $startTimeString);
	fclose($handle);

	fwrite($logOutput, "'sh/hapmap.install_3.php' completed.\n");
	fwrite($logOutput, "...........................................................\n");
	fwrite($logOutput, $currentPath."\n");
	fclose($logOutput);
// Pass control over to a shell script ('sh/hapmap.install_4.sh') to continue processing and link with matlab.
//	$system_call_string = "sh ../sh/hapmap.install_4.sh ".$user." ".$hapmap." > /dev/null &";
	$system_call_string = "sh ../sh/hapmap.install_4.sh ".$user." ".$referencePloidy." ".$project1." ".$project2." ".$hapmap." > /dev/null &";
	system($system_call_string);
?>
<script type="text/javascript">
	var el3 = parent.document.getElementById('Hidden_GenerateNewHapmap');
	el3.style.display = 'none';

	var ff = parent.document.getElementById('panel_hapmap_iframe');
	ff.src = ff.src;
</script>
