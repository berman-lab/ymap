<?php
	session_start();
	error_reporting(E_ALL);
        require_once '../../constants.php';
        ini_set('display_errors', 1);

	if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script> <?php }
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";
	error_reporting(E_ALL);
	ini_set('display_errors', 1);

	// define some characters which shouldn't be in strings.
	$bad_chars          = array(".", ",", "\\", "/", " ", "'", '"', "(", ")", "[", "]");
	$bad_chars2         = array(".", ",", "\\", "/", "'", '"',);

	// Process data transfered to this page.
	$user               = $_SESSION['user'];
	$hapmap             = str_replace($bad_chars, "_",trim( filter_input(INPUT_POST, "hapmap",             FILTER_SANITIZE_STRING) ));
	$genome             = str_replace($bad_chars, "_",trim( filter_input(INPUT_POST, "genome",             FILTER_SANITIZE_STRING) ));
	$project1           = str_replace($bad_chars, "_",trim( filter_input(INPUT_POST, "project1",           FILTER_SANITIZE_STRING) ));
	$project2           = str_replace($bad_chars, "_",trim( filter_input(INPUT_POST, "project2",           FILTER_SANITIZE_STRING) ));
	$colorA             = str_replace($bad_chars, "_",trim( filter_input(INPUT_POST, "homolog_a_color",    FILTER_SANITIZE_STRING) ));
	$colorB             = str_replace($bad_chars, "_",trim( filter_input(INPUT_POST, "homolog_b_color",    FILTER_SANITIZE_STRING) ));
	$hapmap_description = str_replace($bad_chars2,"_",trim( filter_input(INPUT_POST, "hapmap_description", FILTER_SANITIZE_STRING) ));

	// Directories to use later.
	$dir1        = "../../users/".$user."/hapmaps";
	$dir2        = "../../users/".$user."/hapmaps/".$hapmap;
	$dir3        = "../../users/default/hapmaps/".$hapmap;
	$currentPath = getcwd();

	// Initialize 'process_log.txt' file.
	$logOutputName = "../../users/".$user."/hapmaps/".$hapmap."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Running 'scripts_seqModules/scripts_hapmaps/hapmap.update_3.php'.\n");
	fwrite($logOutput, "\tuser            = ".$user."\n");
	fwrite($logOutput, "\thapmap          = ".$hapmap."\n");
	fwrite($logOutput, "\tgenome          = ".$genome."\n");
	fwrite($logOutput, "\tproject1        = ".$project1."\n");
	fwrite($logOutput, "\tproject2        = ".$project2."\n");
	fwrite($logOutput, "\tcolorA          = ".$colorA."\n");
	fwrite($logOutput, "\tcolorB          = ".$colorB."\n");

	// Open 'haplotype.txt' file holding haplotype entry descriptions.
	// Add new haplotype entry descriptions.
	$haplotypeFileName = "../../users/".$user."/hapmaps/".$hapmap."/haplotypeMap.txt";
	$haplotypeFile     = fopen($haplotypeFileName, 'a');
	fwrite($haplotypeFile, $project1."\n".$project2."\n".$hapmap_description."\n");
	fclose($haplotypeFile);

	// Generate 'working.txt' file to let pipeline know that processing is underway.
	$handleName      = "../../users/".$user."/hapmaps/".$hapmap."/working.txt";
	$handle          = fopen($handleName, 'w');
	$startTimeString = date("Y-m-d H:i:s");
	fwrite($handle, $startTimeString);
	fclose($handle);

	fwrite($logOutput, "'scripts_seqModules/scripts_hapmaps/hapmap.update_3.php' completed.\n");
	fwrite($logOutput, "...........................................................\n");
	fwrite($logOutput, $currentPath."\n");
	fclose($logOutput);

	// Pass control over to a shell script ('scripts_seqModules/scripts_hapmaps/hapmap.update_4.sh') to continue processing and link with matlab.
	$system_call_string = "sh hapmap.update_4.sh ".$user." ".$project1." ".$project2." ".$hapmap." > /dev/null &";

	system($system_call_string);
?>
<script type="text/javascript">
	// Close user interface window involved in updating hapmap.
	var el3 = parent.document.getElementById('Hidden_AddToHapmap');
	el3.style.display = 'none';

	// Force reload of hapmap tab in main interface?
	var ff = parent.document.getElementById('panel_hapmap_iframe');
	ff.src = ff.src;
</script>
