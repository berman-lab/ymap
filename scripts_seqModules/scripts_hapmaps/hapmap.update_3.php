<?php
	session_start();
	error_reporting(E_ALL);
        require_once '../../constants.php';
        ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		?> <script type="text/javascript"> parent.reload(); </script> <?php
	}

	// Load user string from session.
	$user       = $_SESSION['user'];

	// Validate hapmap input string.
	$hapmap     = trim(filter_input(INPUT_POST, "hapmap", FILTER_SANITIZE_STRING)); // strip out any html tags.
	$hapmap     = str_replace(" ","_",$hapmap);                                     // convert any spaces to underlines.
	$hapmap     = preg_replace("/[\s\W]+/", "", $hapmap);                           // remove everything but alphanumeric characters and underlines.
	// Confirm if requested hapmap exists.
	$hapmap_dir = "../../users/".$user."/hapmaps/".$hapmap;
	if !(is_dir($hapmap_dir)) {
		// Hapmap doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	// Validate genome input string.
	$genome     = trim(filter_input(INPUT_POST, "genome", FILTER_SANITIZE_STRING));
	$genome     = str_replace(" ","_",$genome);
	$genome     = preg_replace("/[\s\W]+/", "", $genome);
	// Confirm if requested genome exists.
	$genome_dir1 = "../../users/".$user."/genomes/".$genome;
	$genome_dir2 = "../../users/default/genomes/".$genome;
	if !(is_dir($genome_dir1) || is_dir($genome_dir2)) {
		// Genome doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	// Validate project1 input string.
	$project1      = trim(filter_input(INPUT_POST, "project1", FILTER_SANITIZE_STRING));
	$project1      = str_replace(" ","_",$project1);
	$project1      = preg_replace("/[\s\W]+/", "", $project1);
	// Confirm if requested project1 project exists.
	$project1_dir1 = "../../users/".$user."/projects/".$project1;
	$project1_dir2 = "../../users/default/projects/".$project1;
	if !(is_dir($project1_dir1) || is_dir($project1_dir2)) {
		// Parent project doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	// Validate project1 input string.
	$project2      = trim(filter_input(INPUT_POST, "project2", FILTER_SANITIZE_STRING));
	$project2      = str_replace(" ","_",$project2);
	$project2      = preg_replace("/[\s\W]+/", "", $project2);
	// Confirm if requested project2 project exists.
	$project2_dir1 = "../../users/".$user."/projects/".$project2;
	$project2_dir2 = "../../users/default/projects/".$project2;
	if !(is_dir($project2_dir1) || is_dir($project2_dir2)) {
		// Project2 project doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}



	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";

	// Process data transfered to this page.
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
