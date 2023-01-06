<?php
	session_start();
	$user = $_SESSION['user'];

	// This script is intended to take information from file uploaders and then initiate the pipeline scripts to start processing.
	// This has been added to remove any server-side location information from being passed through client-side scripting.


$bad_chars  = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/","'",'"',"[","]","!");
$dataFormat = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "dataFormat", FILTER_SANITIZE_STRING)));
$genome     = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "genome", FILTER_SANITIZE_STRING)));
$project    = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "project", FILTER_SANITIZE_STRING)));
$key        = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "key", FILTER_SANITIZE_STRING)));

// filenames can have characters that shouldn't be in above variables.
$fileName   = filter_input(INPUT_POST, "fileName",   FILTER_SANITIZE_STRING);

$_SESSION['dataFormat'] = $dataFormat;
$_SESSION['fileName']   = $fileName;
$_SESSION['genome']     = $genome;
$_SESSION['project']    = $project;
$_SESSION['key']        = $key;

if (project != "") {
	// initiate project.
	$conclusion_script = "";
	switch ($dataFormat) {
		case "SnpCghArray":
			$conclusion_script = "scripts_SnpCghArray/project.SnpCgh.install.php";
			break;
		case "WGseq_single":
			$conclusion_script = "scripts_seqModules/scripts_WGseq/project.single_WGseq.install_1.php";
			break;
		case "WGseq_paired":
			$conclusion_script = "scripts_seqModules/scripts_WGseq/project.paired_WGseq.install_1.php";
			break;
		case "ddRADseq_single":
			$conclusion_script = "scripts_seqModules/scripts_ddRADseq/project.single_ddRADseq.install_1.php";
			break;
		case "ddRADseq_paired":
			$conclusion_script = "scripts_seqModules/scripts_ddRADseq/project.paired_ddRADseq.install_1.php";
			break;
		case "RNAseq_single":
			$conclusion_script = "scripts_seqModules/scripts_RNAseq/project.single_RNAseq.install_1.php";
			break;
		case "RNAseq_paired":
			$conclusion_script = "scripts_seqModules/scripts_RNAseq/project.paired_RNAseq.install_1.php";
			break;
		case "IonExpressSeq_single":
			$conclusion_script = "scripts_seqModules/scripts_IonExpressSeq/project.single_IonExpressSeq.install_1.php";
			break;
		case "IonExpressSeq_paired":
			$conclusion_script = "scripts_seqModules/scripts_IonExpressSeq/project.paired_IonExpressSeq.install_1.php";
			break;
	}
} else if (genome != "") {
	// initiate genome.
	$conclusion_script = "scripts_genomes/genome.install_1.php";
}
// troubleshooting output
//print "[upload_processer.php]\n";       print "user:        ".$user."\n";              print "project:     ".$project."\n";  print "genome:      ".$genome."\n";
//print "data format: ".$dataFormat."\n"; print "script:      ".$conclusion_script."\n"; print "filename:    ".$fileName."\n"; print "key:         ".$key."\n";

// Move to user directory
chdir("users/".$user);

// Open processing script.
header("Location: ".$conclusion_script);
?>
