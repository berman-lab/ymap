<?php
	session_start();
	error_reporting(E_ALL);
        require_once '../constants.php';
        ini_set('display_errors', 1);

        // If the user is not logged on, redirect to login page.
        if(!isset($_SESSION['logged_on'])){
		session_destroy();
                header('Location: user.login.php');
        }

	// Load user string from session.
	$user   = $_SESSION['user'];

	// Sanitize input strings.
        $key    = trim(filter_input(INPUT_POST, "key", FILTER_SANITIZE_STRING));	// strip out any html tags.
        $key    = str_replace(" ","_",$key);						// convert any spaces to underlines.
        $key    = preg_replace("/[\s\W]+/", "", $key);					// remove everything but alphanumeric characters and underlines.

	// Load strings from session.
	$genome     = $_SESSION['genome_'.$key];
	$chr_count  = $_SESSION['chr_count_'.$key];

	$genome_dir = "../users/".$user."/genomes/".$genome;

?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd"
<HTML>
<HEAD>
	<style type="text/css">
		body {font-family: arial;}
	</style>
<meta http-equiv="content-type" content="text/html; charset=iso-8859-1">
<title>Install genome into pipeline.</title>
</HEAD>
<?php
	$chr_names        = json_decode(file_get_contents($genome_dir."/chr_names.json"), true);
	$chr_lengths      = json_decode(file_get_contents($genome_dir."/chr_lengths.json"), true);;
	$chr_draws        = array();
	$chr_shortNames   = array();
	$chr_cenStarts    = array();
	$chr_cenEnds      = array();
	$chr_count_used   = 0;

// Open 'process_log.txt' file.
	$logOutputName = $genome_dir."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Running 'scripts_genomes/genome.install_2.php'.\n");

	$sizeFile_1   = $genome_dir."/upload_size_1.txt";
	$handle       = fopen($sizeFile_1,'r');
	$sizeString_1 = trim(fgets($handle));
	fclose($handle);
	if ($sizeString_1 !== "") {
		echo "\n<script type='text/javascript'>\n";
		echo "parent.parent.update_genome_file_size('".$key."','".$sizeString_1."');";
		echo "\n</script>\n";
	}

// process POST data.
	fwrite($logOutput, "\tProcessing POST data containing genome specific information.\n");
	$rDNA_start         = filter_input(INPUT_POST, "rDNAstart",          FILTER_SANITIZE_STRING);
	$rDNA_end           = filter_input(INPUT_POST, "rDNAend",            FILTER_SANITIZE_STRING);
	$ploidyDefault      = filter_input(INPUT_POST, "ploidy",             FILTER_SANITIZE_STRING);
	$annotation_count   = filter_input(INPUT_POST, "annotation_count",   FILTER_SANITIZE_STRING);
	$expression_regions = filter_input(INPUT_POST, "expression_regions", FILTER_SANITIZE_STRING);

	if ($chr_count != 0) {
		for ($chr=0; $chr<$chr_count; $chr += 1) {
			$chrID = $chr + 1;
			if (filter_input(INPUT_POST, "draw_".$chrID, FILTER_SANITIZE_STRING) == "on") {
				$chr_draw       = 1;
				$chr_count_used += 1;
			} else {
				$chr_draw = 0;
			}
			$chr_shortName        = filter_input(INPUT_POST, "short_".$chrID,    FILTER_SANITIZE_STRING);
			$chr_cenStart         = filter_input(INPUT_POST, "cenStart_".$chrID, FILTER_SANITIZE_STRING);
			$chr_cenEnd           = filter_input(INPUT_POST, "cenEnd_".$chrID,   FILTER_SANITIZE_STRING);
			$chr_draws[$chr]      = $chr_draw;
			$chr_shortNames[$chr] = $chr_shortName;
			$chr_cenStarts[$chr]  = $chr_cenStart;
			$chr_cenEnds[$chr]    = $chr_cenEnd;
		}
	}
	if (isset($_POST['rDNAchr']) && !empty($_POST['rDNAchr'])) {
		$rDNA_chr = filter_input(INPUT_POST, "rDNAchr", FILTER_SANITIZE_STRING);
	} else {
		$rDNA_chr = "null";
	}

// Generate 'chromosome_sizes.txt' :
	fwrite($logOutput, "\tGenerating 'chromosome_sizes.txt' file.\n");
	$outputName       = $genome_dir."/chromosome_sizes.txt";
	if (file_exists($outputName)) {
		$fileContents = file_get_contents($outputName);
		unlink($outputName);
		$output       = fopen($outputName, 'w');
		fwrite($output, $fileContents);
	} else {
		$output       = fopen($outputName, 'w');
		fwrite($output, "# Chr\tsize(bp)\tname\n");
		for ($chr=0; $chr<$chr_count; $chr += 1) {
			$chrID = $chr + 1;
			if ($chr_draws[$chr] == 1) {
				fwrite($output, $chrID."\t".$chr_lengths[$chr]."\t".$chr_shortNames[$chr]."\n");
			}
		}
		fclose($output);
	}

// Generate 'centromere_locations.txt' :
	fwrite($logOutput, "\tGenerating 'centromere_locations.txt' file.\n");
	$outputName       = $genome_dir."/centromere_locations.txt";
	if (file_exists($outputName)) {
		$fileContents = file_get_contents($outputName);
		unlink($outputName);
		$output       = fopen($outputName, 'w');
		fwrite($output, $fileContents);
	} else {
		$output       = fopen($outputName, 'w');
		fwrite($output, "# Chr\tCEN-start\tCEN-end\n");
		for ($chr=0; $chr<$chr_count; $chr += 1) {
			$chrID    = $chr + 1;
			if ($chr_draws[$chr] == 1) {
				fwrite($output, $chrID."\t".$chr_cenStarts[$chr]."\t".$chr_cenEnds[$chr]."\n");
			}
		}
	}
	fclose($output);

// Generate 'figure_definitions.txt' :
	fwrite($logOutput, "\tGenerating 'figure_definitions.txt' file.\n");
	$max_length       = max($chr_lengths);
	$outputName       = $genome_dir."/figure_definitions.txt";
	if (file_exists($outputName)) {
		$fileContents = file_get_contents($outputName);
		unlink($outputName);
		$output       = fopen($outputName, 'w');
		fwrite($output, $fileContents);
	} else {
		$output       = fopen($outputName, 'w');
		fwrite($output, "# Chr\tUse\tLabel\tName\tposX\tposY\twidth\theight\n");
		if ($chr_count != 0) {
			$usedChrID = 0; // used to count the number of used chromosomes that will be drawn for positioning of the stacked figure
			// setting figure height to be the same for all figures making them ocuppy 50 precent of the maximum height (50 precent for gap)
			$fig_height = 0.5*(0.97/($chr_count_used + 0.5));
			for ($chr=0; $chr<$chr_count; $chr += 1) {
				$chrID = $chr + 1;
				if ($chr_draws[$chr] == 1) {// if this chromosome should be drawn incrementing
				  $usedChrID += 1;
				}
				$fig_posX   = 0.15;
				// title gets 0.03 of the space, and figures share the rest (+0.5 to avoid cutting in the end)
				$fig_posY   = 0.97-(0.97/($chr_count_used + 0.5))*$usedChrID;
				if ($chr_lengths[$chr] == $max_length) {
					$fig_width = "0.8";
				} else {
					$fig_width = "*";
				}
				if ($chr_draws[$chr] == 1) {
					fwrite($output, $chrID."\t1\t".$chr_shortNames[$chr]."\t".$chr_names[$chr]."\t".$fig_posX."\t".$fig_posY."\t".$fig_width."\t".$fig_height."\n");
				} else {
					fwrite($output, "0\t0\t".$chr_shortNames[$chr]."\t".$chr_names[$chr]."\t0\t0\t0\t0\n");
				}
			}
		}
	}
	fclose($output);

// Generate "genome.bed" :
	fwrite($logOutput, "\tGenerating 'genome.bed' file.\n");
	$max_length       = max($chr_lengths);
	$outputName       = $genome_dir."/genome.bed";
	if (file_exists($outputName)) {
		$fileContents = file_get_contents($outputName);
		unlink($outputName);
		$output       = fopen($outputName, 'w');
		fwrite($output, $fileContents);
	} else {
		$output       = fopen($outputName, 'w');
		fwrite($output, "# BED file to limit Abr2 to only indel-realigning used chromosomes.\n");
		fwrite($output, "# [chr name]\t0\t[chr length]\n");
		if ($chr_count != 0) {
			for ($chr=0; $chr<$chr_count; $chr += 1) {
				$chrID = $chr + 1;
				if ($chr_draws[$chr] == 1) {
					fwrite($output, $chr_names[$chr]."\t0\t".$chr_lengths[$chr]."\n");
				}
			}
		}
	}
	fclose($output);

// Generate "ploidy.txt" :
	fwrite($logOutput, "\tGenerating 'ploidy.txt' file.\n");
	$outputName       = $genome_dir."/ploidy.txt";
	if (file_exists($outputName)) {
		$fileContents = file_get_contents($outputName);
		unlink($outputName);
		$output       = fopen($outputName, 'w');
		fwrite($output, $fileContents);
	} else {
		$output       = fopen($outputName, 'w');
		fwrite($output, $ploidyDefault);
	}
	fclose($output);

	// Save important variables to $_SESSION.
	fwrite($logOutput, "\tStoring PHP session variables.\n");
	$_SESSION['rDNA_chr_'.$key]         = $rDNA_chr;
	$_SESSION['rDNA_start_'.$key]       = $rDNA_start;
	$_SESSION['rDNA_end_'.$key]         = $rDNA_end;
	$_SESSION['ploidyDefault_'.$key]    = $ploidyDefault;
	$_SESSION['annotation_count_'.$key] = $annotation_count;

// Debugging output of all variables.
//	print_r($GLOBALS);
	fwrite($logOutput, "\tGenerating form to request annotation information from the user.\n");
?>
<BODY>
	<font color="red" size="2">Enter genome annotation details:</font>
	<form name="annotation_form" id="annotation_form" action="genome.install_3.php" method="post">
		<table border="0">
		<tr>
			<th rowspan="2"><font size="2">Chr</font></th>
			<th rowspan="2"><font size="2">Type</font></th>
			<th rowspan="2"><font size="2">Start BP</font></th>
			<th rowspan="2"><font size="2">End BP</font></th>
			<th rowspan="2"><font size="2">Name</font></th>
			<th colspan="2"><font size="2">Color</font></th>
			<th rowspan="2"><font size="2">Size</font></th>
		</tr>
		<tr>
			<th><font size="2">Fill</font></th>
			<th><font size="2">Edge</font></th>
		</tr>
<?php	if ($rDNA_chr != "null") {?>
		<tr>
			<td align="middle"><font size="2"><?php echo $chr_shortNames[$rDNA_chr-1]; ?></font></td>
			<td align="middle"><font size="2">dot</font></td>
			<td align="middle"><font size="2"><?php echo $rDNA_start; ?></font></td>
			<td align="middle"><font size="2"><?php echo $rDNA_end; ?></font></td>
			<td align="middle"><font size="2">rDNA</font></td>
			<td align="middle"><font size="2">blue</font></td>
			<td align="middle"><font size="2">blue</font></td>
			<td align="middle"><font size="2">5</font></td>
		</tr>
<?php
		}
			for ($annotation=0; $annotation<$annotation_count; $annotation+=1) {
				echo "\t\t<tr>\n";
				echo "\t\t\t<td><select name=\"annotation_chr_{$annotation}\">";
				for ($chr=0; $chr<$chr_count; $chr += 1) {
					$chrID = $chr + 1;
					if ($chr_draws[$chr] == 1) {
						echo "\t\t\t\t<option value=\"{$chrID}\">".$chr_shortNames[$chr]."</option>";
					}
				}
				echo "</select>\t\t\t</td>\n";

				echo "\t\t\t<td><select name=\"annotation_shape_{$annotation}\">";
				echo "\t\t\t\t<option value=\"dot\">dot</option>";
				echo "\t\t\t\t<option value=\"block\">block</option>";
				echo "</select>\t\t\t</td>\n";

				echo "\t\t\t<td><input type=\"text\" name=\"annotation_start_{$annotation}\" value=\"0\" size=\"6\"></td>\n";
				echo "\t\t\t<td><input type=\"text\" name=\"annotation_end_{$annotation}\" value=\"0\" size=\"6\"></td>\n";
				echo "\t\t\t<td><input type=\"text\" name=\"annotation_name_{$annotation}\" value=\"0\" size=\"6\"></td>\n";

				echo "\t\t\t<td><select name=\"annotation_fillColor_{$annotation}\">";
				echo "\t\t\t\t<option value=\"k\">black</option>";
				echo "\t\t\t\t<option value=\"y\">yellow</option>";
				echo "\t\t\t\t<option value=\"m\">magenta</option>";
				echo "\t\t\t\t<option value=\"c\">cyan</option>";
				echo "\t\t\t\t<option value=\"r\">red</option>";
				echo "\t\t\t\t<option value=\"g\">green</option>";
				echo "\t\t\t\t<option value=\"b\">blue</option>";
				echo "\t\t\t\t<option value=\"w\">white</option>";
				echo "</select>\t\t\t</td>\n";

				echo "\t\t\t<td><select name=\"annotation_edgeColor_{$annotation}\">";
				echo "\t\t\t\t<option value=\"k\">black</option>";
				echo "\t\t\t\t<option value=\"y\">yellow</option>";
				echo "\t\t\t\t<option value=\"m\">magenta</option>";
				echo "\t\t\t\t<option value=\"c\">cyan</option>";
				echo "\t\t\t\t<option value=\"r\">red</option>";
				echo "\t\t\t\t<option value=\"g\">green</option>";
				echo "\t\t\t\t<option value=\"b\">blue</option>";
				echo "\t\t\t\t<option value=\"w\">white</option>";
				echo "</select>\t\t\t</td>\n";

				echo "\t\t\t<td><input type=\"text\" name=\"annotation_size_{$annotation}\" value=\"5\" size=\"6\"></td>\n";

				echo "\t\t</tr>\n";
			}
		?>
		</table><br>
		<input type="submit" id="form_submit_gi2" value="Save genome details..." disabled>

		<input type="hidden" id="key"                name="key"                value="<?php echo $key;                ?>">
		<input type="hidden" id="expression_regions" name="expression_regions" value="<?php echo $expression_regions; ?>">
	</form>
<!------------------ Onload-script section ------------------!>
	<script type="text/javascript">
	onload=function(){
		if ("<?php echo $annotation_count; ?>" == "") {
			document.getElementById('form_submit_gi2').disabled = false;
			document.annotation_form.submit();
		} else {
			if (<?php echo $annotation_count; ?> == 0) {
				document.getElementById('form_submit_gi2').disabled = false;
				document.annotation_form.submit();
			} else {
				parent.parent.resize_genome('<?php echo $key; ?>', <?php echo ($annotation_count + 5.5)*27; ?>);
				document.getElementById('form_submit_gi2').disabled = false;
			}
		}
	}
	</script>
</BODY>
</HTML>
<?php
	fwrite($logOutput, "\t'scripts_genomes/genome_install_2.php' has completed.\n");
	fclose($logOutput);
?>
