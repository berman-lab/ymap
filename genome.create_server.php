<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		header('Location: user.login.php');
	}

//	// If the user is not logged on, redirect to login page.
//	if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script> <?php }

//      $user          = $_SESSION['user'];
//      $genome        = trim(filter_input(INPUT_POST, "genome", FILTER_SANITIZE_STRING));  // removing unwanted characters

	$bad_chars      = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/","'",'"',"[","]","!");
	$genomeName = trim(filter_input(INPUT_POST, "newGenomeName", FILTER_SANITIZE_STRING));	// removing unwanted characters
	$genomeNameTrimmed     = str_replace($bad_chars,"",$genomeName);			// changing spaces to underlines
	$genomeNameTrimmed     = str_replace(" ","_",$genomeNameTrimmed);
	$user           = $_SESSION['user'];
	$dir1           = "users/".$user."/genomes";
	$dir2           = "users/".$user."/genomes/".$genomeNameTrimmed;
	$dir3           = "users/default/genomes/".$genomeNameTrimmed;


	// Deals with accidental deletion of genomes dir.
	if (!file_exists($dir1)){
		mkdir($dir1);
		chmod($dir1,0777);
	}

	if (file_exists($dir2) || file_exists($dir3)) {
		// Directory already exists
		echo "Genome '".$genomeName."' directory already exists.";
?>
	<html>
	<body>
	<script type="text/javascript">
	var el1 = parent.document.getElementById('Hidden_InstallNewGenome');
	el1.style.display = 'none';

	var el2 = parent.document.getElementById('panel_genome_iframe').contentDocument.getElementById('name_error_comment');
	el2.style.visibility = 'visible';

	window.location = "genome.create_window.php";
	</script>
	</body>
	</html>
<?php
	} else {
		// Create the genome folder inside the user's genomes directory
		mkdir($dir2);
		chmod($dir2,0777);

		// Generate 'name.txt' file containing:
		//      one line; name of genome.
		$outputName   = "users/".$user."/genomes/".$genomeNameTrimmed."/name.txt";
		$output       = fopen($outputName, 'w');
		fwrite($output, $genomeName);
		fclose($output);

		$_SESSION['pending_install_genome_count'] += 1;
?>
	<html>
	<body>
	<script type="text/javascript">
	var el1 = parent.document.getElementById('panel_genome_iframe').contentDocument.getElementById('newly_installed_list');
	el1.innerHTML += "<?php echo $_SESSION['pending_install_genome_count']; ?>. <?php echo $genomeName; ?><br>";

	var el2 = parent.document.getElementById('panel_genome_iframe').contentDocument.getElementById('pending_comment');
	el2.style.visibility = 'visible';

	var el3 = parent.document.getElementById('panel_genome_iframe').contentDocument.getElementById('name_error_comment');
	el3.style.visibility = 'hidden';

	var el4 = parent.document.getElementById('Hidden_InstallNewGenome');
	el4.style.display = 'none';

	window.location = "genome.create_window.php";
	</script>
	</body>
	</html>
<?php
	}
?>
