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
	$user   = $_SESSION['user'];

	// Validate hapmap input string.
	$hapmap = trim(filter_input(INPUT_POST, "hapmap", FILTER_SANITIZE_STRING));	// strip out any html tags.
	$hapmap = str_replace(" ","_",$hapmap);						// convert any spaces to underlines.
	$hapmap = preg_replace("/[\s\W]+/", "", $hapmap);				// remove everything but alphanumeric characters and underlines.
	// Confirm if requested genome exists.
	$hapmap_dir         = "../../users/".$user."/hapmaps/".$hapmap;
	$default_hapmap_dir = "../../users/default/hapmaps/".$hapmap;
	if (!is_dir($hapmap_dir)) {
		// Hapmap doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	// Validate genome input string.
	$genome = trim(filter_input(INPUT_POST, "genome", FILTER_SANITIZE_STRING));
	$genome = str_replace(" ","_",$genome);
	$genome = preg_replace("/[\s\W]+/", "", $genome);
	// Confirm if requested genome exists.
	$genome_dir = "../../users/".$user."/genomes/".$genome;
	if (!is_dir($genome_dir)) {
		// Genome doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	// Validate project input string.
	$project1        = trim(filter_input(INPUT_POST, "project1", FILTER_SANITIZE_STRING));	// strip out any html tags.
	$project1        = str_replace(" ","_",$project1);					// convert any spaces to underlines.
	$project1        = preg_replace("/[\s\W]+/", "", $project1);				// remove everything but alphanumeric characters and underlines.
	$project2        = trim(filter_input(INPUT_POST, "project2", FILTER_SANITIZE_STRING));
	$project2        = str_replace(" ","_",$project2);
	$project2        = preg_replace("/[\s\W]+/", "", $project2);
	// Confirm if requested projects exists.
	$project1_dir1   = "users/".$user."/projects/".$project1;
	$project1_dir2   = "users/default/projects/".$project1;
	$project2_dir1   = "users/".$user."/projects/".$project2;
	$project2_dir2   = "users/default/projects/".$project2;
	if !(is_dir($project1_dir1) || is_dir($project1_dir2) || is_dir($project2_dir1) || is_dir($project2_dir2)) {
		// A project doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	// Validate color input strings
	$colorA          = trim(filter_input(INPUT_POST, "homolog_a_color", FILTER_SANITIZE_STRING));
	$colorA          = str_replace(" ","_",$colorA);
	$colorA          = preg_replace("/[\s\W]+/", "", $colorA);
	$colorB          = trim(filter_input(INPUT_POST, "homolog_b_color", FILTER_SANITIZE_STRING));
	$colorB          = str_replace(" ","_",$colorB);
	$colorB          = preg_replace("/[\s\W]+/", "", $colorB);

	// Validate reference ploidy input string.
        $referencePloidy = trim(filter_input(INPUT_POST, "referencePloidy", FILTER_SANITIZE_STRING));	// strip out any html tags.
	$referencePloidy = (float)preg_replace("/[^\d\.]+/", "", $referencePloidy);			// remove everything but numerals and period.
	if ($referencePloidy == 2) {
		// Validate hapmap description string.
		$hapmap_description = trim(filter_input(INPUT_POST, "hapmap_description", FILTER_SANITIZE_STRING));
	}

	$currentPath = getcwd();

        if (file_exists($project_dir) || file_exists($default_project_dir)) {
		//============================================
		// Hapmap directory already exists, so exit.
		//--------------------------------------------
		echo "Hapmap '".$hapmap."' directory already exists.";
?>
        <html>
        <body>
        <script type="text/javascript">
        var el1 = parent.document.getElementById('Hidden_InstallNewDataset');
        el1.style.display = 'none';

        var el2 = parent.document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('name_error_comment');
        el2.style.visibility = 'visible';

        window.location = "project.create_window.php";
        </script>
        </body>
        </html>
<?php
	} else {
		//========================================================
		// Hapmap directory doesn't exist, go about creating it.
		//--------------------------------------------------------

		// Create the hapmap folder inside the user's hapmaps directory
		mkdir($hapmap_dir);
		chmod($hapmap_dir,0777);

		// Initialize 'process_log.txt' file.
		$logOutputName = $hapmap_dir."/process_log.txt";
		$logOutput     = fopen($logOutputName, 'a');
		fwrite($logOutput, "Running 'scripts_seqModules/scripts_hapmaps/hapmap.install_3.php'.\n");
		fwrite($logOutput, "\tuser            = ".$user."\n");
		fwrite($logOutput, "\thapmap          = ".$hapmap."\n");
		fwrite($logOutput, "\tgenome          = ".$genome."\n");
		fwrite($logOutput, "\treferencePloidy = ".$referencePloidy."\n");
		fwrite($logOutput, "\tproject1        = ".$project1."\n");
		fwrite($logOutput, "\tproject2        = ".$project2."\n");
		fwrite($logOutput, "\tcolorA          = ".$colorA."\n");
		fwrite($logOutput, "\tcolorB          = ".$colorB."\n");

		// Create 'colors.txt' file to contain colors used in haplotype figures.
		$handleName = $hapmap_dir."/colors.txt";
		if (file_exists($handleName)) {
		} else {
			$handle     = fopen($handleName, 'w');
			fwrite($handle, $colorA."\n".$colorB);
			fclose($handle);
		}

		// Create 'genome.txt' file to contain genome used in haplotype.
		$handleName = $hapmap_dir."/genome.txt";
		if (file_exists($handleName)) {
		} else {
			$handle     = fopen($handleName, 'w');
			fwrite($handle, $genome);
			fclose($handle);
		}

		// Create 'parent.txt' file to contain parent genome used in haplotype.
		$handleName = $hapmap_dir."/parent.txt";
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
			$haplotypeFileName = $hapmap_dir."/haplotypeMap.txt";
			$haplotypeFile     = fopen($haplotypeFileName, 'a');
			fwrite($haplotypeFile, $project1."\n".$project2."\n".$hapmap_description."\n");
			fclose($haplotypeFile);
		} else {
			// for haplotype map derived from two haploid strains, there will not be more entries after construction.
			// Because of this, 'parent.txt' will contain sufficient information to define haplotype map.
		}

		// Generate 'working.txt' file to let pipeline know that processing is underway.
		$handleName      = $hapmap_dir."/working.txt";
		$handle          = fopen($handleName, 'w');
		$startTimeString = date("Y-m-d H:i:s");
		fwrite($handle, $startTimeString);
		fclose($handle);

		fwrite($logOutput, "'scripts_seqModules/scripts_hapmaps/hapmap.install_3.php' completed.\n");
		fwrite($logOutput, "...........................................................\n");
		fwrite($logOutput, $currentPath."\n");
		fclose($logOutput);

		// Pass control over to a shell script ('scripts_seqModules/scripts_hapmaps/hapmap.install_4.sh') to continue processing and link with matlab.
		$system_call_string = "sh hapmap.install_4.sh ".$user." ".$referencePloidy." ".$project1." ".$project2." ".$hapmap." > /dev/null &";
		system($system_call_string);
?>
<script type="text/javascript">
	var el3 = parent.document.getElementById('Hidden_GenerateNewHapmap');
	el3.style.display = 'none';

	var ff = parent.document.getElementById('panel_hapmap_iframe');
	ff.src = ff.src;
</script>
<?php
        }
?>
