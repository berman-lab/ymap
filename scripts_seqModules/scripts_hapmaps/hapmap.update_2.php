<?php
	session_start();
	if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script> <?php } else { $user = $_SESSION['user']; }
	require_once '../../constants.php';
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";

	error_reporting(E_ALL);
	ini_set('display_errors', 1);

	// define some characters which shouldn't be in strings.
	$bad_chars       = array(".", ",", "\\", "/", " ", "'", '"', "(", ")", "[", "]");

	$user            = $_SESSION['user'];
	$hapmap          = str_replace($bad_chars,"_",trim( filter_input(INPUT_POST, "hapmap",   FILTER_SANITIZE_STRING) ));
	$genome          = str_replace($bad_chars,"_",trim( filter_input(INPUT_POST, "genome",   FILTER_SANITIZE_STRING) ));
	$project1        = str_replace($bad_chars,"_",trim( filter_input(INPUT_POST, "project1", FILTER_SANITIZE_STRING) ));
	$project2        = str_replace($bad_chars,"_",trim( filter_input(INPUT_POST, "project2", FILTER_SANITIZE_STRING) ));

	// Load "chromosome_sizes.txt" file from in-use genome.
	$file1 = "../../users/default/genomes/".$genome."/chromosome_sizes.txt";
	$file2 = "../../users/".$user."/genomes/".$genome."/chromosome_sizes.txt";
	if (file_exists($file1)) {
		$handle     = fopen($file1,'r');
		$fileString = fread($handle, filesize($file1));
	} else if (file_exists($file2)) {
		$handle     = fopen($file2,'r');
		$fileString = fread($handle, filesize($file2));
	} else {
		exit;
	}
	fclose($handle);

	// Determine number of chromosomes.
	$chrLinesArray  = explode("\n", trim($fileString));
	array_shift($chrLinesArray);
	$chr_count      = count($chrLinesArray);

	// Build array of chromosome information.
	$chrNumID       = array();
	$chrNameShort   = array();
	$chrSize        = array();
	$newTags        = array();
	for ($chr=0; $chr<$chr_count; $chr+=1) { // Build out arrays with empty values.
		array_push($chrNumID,"");
		array_push($chrNameShort,"");
		array_push($chrSize,"");
		array_push($newTags,"");
	}
	foreach ($chrLinesArray as $key=>$chrLine) { // Populate arrays with data.
		$chrLineParts       = preg_split('/\s+/', $chrLine);
		$oneChrNum          = $chrLineParts[0];
		$oneChrSize         = $chrLineParts[1];
		$oneChrName         = $chrLineParts[2];
		$chrNumID[$key]     = $oneChrNum;
		$chrNameShort[$key] = $oneChrName;
		$chrSize[$key]      = $oneChrSize;
	}

	// dragon
	for ($chr=0; $chr<$chr_count; $chr+=1) {
		$newTags[$chr] = filter_input(INPUT_POST, "newTag_{$chr}", FILTER_SANITIZE_STRING);
	}
	$newTags[0] = filter_input(INPUT_POST, "newTag", FILTER_SANITIZE_STRING);
	print_r($newTags[0]);

	echo "<script type='text/javascript'>\n";
	echo "chr_count = ".$chr_count."\n";
	echo "</script>\n";
?>
<HEAD>
	<style type="text/css">
		body {font-family: arial;}
		.tab {margin-left:   1cm;}
		select.list1 option.optionCyan {background-color: #007700;}
	</style>
	<meta http-equiv="content-type" content="text/html; charset=utf-8">
	<title>[Needs Title]</title>
</HEAD>
<BODY onload="showColors('homolog_a_color','colorBar1'); showColors('homolog_b_color','colorBar2')">
<b>Fill in details for this hapmap entry.</b>
<div class="tab">
	Genome : <?php echo $genome; ?><br>
	<?php
	// Find and display parent strain figure.
	echo "Reference dataset : ".$project1."<br>";
	if (file_exists("../../users/".$user."/projects/".$project1)) {
	    if (file_exists("../../users/".$user."/projects/".$project1."/fig.CNV-LOH-map.2.png")) {
	        $imageUrl = "../../users/".$user."/projects/".$project1."/fig.CNV-LOH-map.2.png";
	    } else {
	        $imageUrl = "../../users/".$user."/projects/".$project1."/fig.CNV-SNP-map.2.png";
	    }
	    echo "<img src=\"{$imageUrl}\" width=\"50%\">\n";
	} else if (file_exists("../../users/default/projects/".$project1)) {
		if (file_exists("../../users/default/projects/".$project1."/fig.CNV-LOH-map.2.png")) {
			$imageUrl = "../../users/default/projects/".$project1."/fig.CNV-LOH-map.2.png";
		} else {
			$imageUrl = "../../users/default/projects/".$project1."/fig.CNV-SNP-map.2.png";
		}
		echo "<img src=\"{$imageUrl}\" width=\"50%\">\n";
	} else {
		exit;
	}
	echo "<br>";
	// Find and display secondary strain figure.
	echo "Experimental dataset : ".$project2."<br>";
	if (file_exists("../../users/".$user."/projects/".$project2)) {
		if (file_exists("../../users/".$user."/projects/".$project2."/fig.CNV-LOH-map.2.png")) {
			$imageUrl = "../../users/".$user."/projects/".$project2."/fig.CNV-LOH-map.2.png";
		} else {
			$imageUrl = "../../users/".$user."/projects/".$project2."/fig.CNV-SNP-map.2.png";
		}
		echo "<img src=\"{$imageUrl}\" width=\"50%\">\n";
		$CGD_annotations_url = "../../users/".$user."/projects/".$project2."/CGD_annotations.".$project2.".txt";
		if (file_exists("../../users/".$user."/projects/".$project2."/CGD_annotations.".$project2.".txt")) {
			echo "<br><div class='tab'>Examine <button onclick=\"loadExternal('".$CGD_annotations_url."',50,  220);\">GBrowse annotation track</button> to determine precise breakpoints.</div><br>";
		}
	} else if (file_exists("../../users/default/projects/".$project2)) {
		if (file_exists("../../users/default/projects/".$project2."/fig.CNV-LOH-map.2.png")) {
			$imageUrl = "../../users/default/projects/".$project2."/fig.CNV-LOH-map.2.png";
		} else {
			$imageUrl = "../../users/default/projects/".$project2."/fig.CNV-SNP-map.2.png";
		}
		echo "<img src=\"{$imageUrl}\" width=\"50%\">\n";
		$CGD_annotations_url = "../../users/default/projects/".$project2."/CGD_annotations.".$project2.".txt";
		if (file_exists("../../users/default/projects/".$project2."/CGD_annotations.".$project2.".txt")) {
			echo "<button onclick=\"loadExternal('".$CGD_annotations_url."',50,  220);\">GBrowse</button><br>";
		}
	} else {
		exit;
	}
	?>
	<br>
</div>
<table><tr><td>
	<div class="tab">
		<table border="0">
		<tr>
			<th bgcolor=\"#CCFFCC\"><i>Chr</i></th>
			<th bgcolor=\"#CCCCFF\"><i>Homolog</i></th>
			<th bgcolor=\"#CCFFCC\"><i>Start bp</i></th>
			<th bgcolor=\"#CCCCFF\"><i>End bp</i></th>
			<th bgcolor=\"#CCFFCC\"></th>
			<th bgcolor=\"#CCCCFF\"></th>
			<th bgcolor=\"#CCFFCC\"></th>
		</tr>
		<?php
		// generate line in table for each chromosome.
		for ($chr=0; $chr<$chr_count; $chr+=1) {
			$chrID = $chr+1;
			echo "\t\t<tr>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCFFCC\">{$chrNameShort[$chr]}</td>\n";
			echo "\t\t\t                                         <input type=\"hidden\"   id=\"chr_{$chr}\"     value=\"{$chrNumID[$chr]}\"></input>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCCCFF\"><select id=\"homolog_{$chr}\"><option value=\"a\">a</option><option value=\"b\">b</option></select></td>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCFFCC\"><input type=\"text\"     id=\"start_{$chr}\"   value=\"1\"                 size=\"8\"></td>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCCCFF\"><input type=\"text\"     id=\"end_{$chr}\"     value=\"{$chrSize[$chr]}\"  size=\"8\"></td>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCFFCC\"><button onclick=\"updateTag_{$chr}();\">+</button></td>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCCCFF\"><div id=\"tag_{$chr}\">{$newTags[$chr]}</div></td>\n";
			echo "\t\t\t<td align=\"center\" bgcolor=\"#CCFFCC\"><button onclick=\"clearTag_{$chr}();\">-</button></td>\n";
			echo "\t\t</tr>\n";
		}
		?>
		</table>
		<script type="text/javascript">
		<?php
		// Generate javascript for each chromosome.
		for ($chr=0; $chr<$chr_count; $chr+=1) {
			echo "updateTag_{$chr} = function() {\n";
			echo "\tchr_id        = document.getElementById(\"chr_{$chr}\"    ).value;\n";
			echo "\thomolog_id    = document.getElementById(\"homolog_{$chr}\").value;\n";
			echo "\tstart_bp      = document.getElementById(\"start_{$chr}\"  ).value;\n";
			echo "\tend_bp        = document.getElementById(\"end_{$chr}\"    ).value;\n";
			echo "\tnewString     = String(homolog_id)+'['+String(start_bp)+':'+String(end_bp)+'] '\n";
			echo "\tdocument.getElementById(\"tag_{$chr}\").innerHTML += newString;\n";
			echo "}\n";
			echo "clearTag_{$chr} = function() {\n";
	        echo "\tdocument.getElementById(\"tag_{$chr}\").innerHTML = \"\";\n";
	        echo "}\n";
		}
		?>
		function loadExternal(imageUrl) {
			window.open(imageUrl);
		}
		</script>
	</div>
</td><td valign="top">
    <div class="tab">
        <br>
        If the experimental dataset was analyzed using the reference dataset, any regions of loss of heterozygosity (LOH) would appear in red.<br>
        If the strains are unrelated, it is unlikely this analysis will be helpful.<br><br>
        Enter a description of the homolog identity of the homozygosed regions using the form at left.<br>
        <div class="tab">
            Set the '<i><b>Homolog</b></i>', '<i><b>Start bp</b></i>', and '<i><b>End bp</b></i>' options to describe a LOH region on a chromosome.<br>
            Click the <b>[+]</b> button to accept the LOH description.<br>
            Click the <b>[-]</b> buttons to clear the entry for a chromosome.
        </div>
		<br>
		Clicking the '<b>GBrowse</b>' button will open a custom GBrowse annotation track file for the experimental dataset in a new window.   The
		annotation file can be loaded into GBrowse at <a href="http://Candidagenome.org">Candida Genome Database</a>, for supported genomes.   This
		can be useful in determining the end-coordinates of LOH regions.
    </div>
</td></tr>
<tr><td>
<div class="tab">
	<?php
	// Load the colors defined for this hapmap.
	$colorsFile = "../../users/".$user."/hapmaps/".$hapmap."/colors.txt";
	$handle       = fopen($colorsFile,'r');
	$colorString1 = trim(fgets($handle));
	$colorString2 = trim(fgets($handle));
	?>
	<br>
	Homolog 'a' color : <input type="text" id="homolog_a_color" name="homolog_a_color" value="<?php echo $colorString1; ?>" readonly>
	<div id="colorBar1"></div>
	<br>
	Homolog 'b' color : <input type="text" id="homolog_b_color" name="homolog_b_color" value="<?php echo $colorString2; ?>" readonly>
	<div id="colorBar2"></div>
	<br><br>
<div>
<script type="text/javascript">
	// Show the user the colors that are defined for this hapmap.
	// Script run on page load.
	showColors=function(listToLookAt, targetToChange) {
		var selectedColor = document.getElementById(listToLookAt).value;
		var select        = document.getElementById(targetToChange);
		select.innerHTML  = '';
		if        (selectedColor == "deep pink")       { select.innerHTML = '<div class="tab" style="background-color:rgb(255,0  ,127)"> &nbsp;</div>';
		} else if (selectedColor == "magenta")         { select.innerHTML = '<div class="tab" style="background-color:rgb(255,0  ,255)"> &nbsp;</div>';
		} else if (selectedColor == "electric indigo") { select.innerHTML = '<div class="tab" style="background-color:rgb(127,0  ,255)"> &nbsp;</div>';
		} else if (selectedColor == "blue")            { select.innerHTML = '<div class="tab" style="background-color:rgb(0  ,0  ,255)"> &nbsp;</div>';
		} else if (selectedColor == "dodger blue")     { select.innerHTML = '<div class="tab" style="background-color:rgb(0  ,127,255)"> &nbsp;</div>';
		} else if (selectedColor == "cyan")            { select.innerHTML = '<div class="tab" style="background-color:rgb(0  ,255,255)"> &nbsp;</div>';
		} else if (selectedColor == "spring green")    { select.innerHTML = '<div class="tab" style="background-color:rgb(0  ,255,127)"> &nbsp;</div>';
		} else if (selectedColor == "green")           { select.innerHTML = '<div class="tab" style="background-color:rgb(0  ,255,0  )"> &nbsp;</div>';
		} else if (selectedColor == "chartreuse")      { select.innerHTML = '<div class="tab" style="background-color:rgb(127,255,0  )"> &nbsp;</div>';
		} else if (selectedColor == "yellow")          { select.innerHTML = '<div class="tab" style="background-color:rgb(255,255,0  )"> &nbsp;</div>';
		} else if (selectedColor == "dark orange")     { select.innerHTML = '<div class="tab" style="background-color:rgb(255,127,0  )"> &nbsp;</div>';
		} else {                                         select.innerHTML = '<div class="tab" style="background-color:rgb(127,127,127)"> &nbsp;</div>';
		}
	}
</script>

</td><td valign="top">
<?php
	if (file_exists($colorsFile)) {
	} else {
		?>
		<div class="tab">
			<br>
			Select which colors to be used in displaying the homologs for this haplotype map.<br>
			Colors should be chosen to minimize visual confusion.<br>
			<br>
			Red is not available, as it is the color used to highlight LOHs in regions without haplotype definitions.
		</div>
		<?php
	}
?>
</td></tr>
<tr><td>
<div class="tab">
	<button onclick="GatherAndSubmitData();">Process haplotype entry...</button>
	<script type="text/javascript">
		GatherAndSubmitData=function() {
			var string_allData = "";
			for (var chr=0;chr<chr_count; chr++) {
				var string_start   = String(chr)+"(";
				var string_start   = document.getElementById("chr_"+String(chr)).value+"(";
				var string_inner   = trim(document.getElementById("tag_"+String(chr)).innerHTML);
				var string_end     = ")";
				if (string_inner.length > 0) {
					string_allData += string_start+string_inner+string_end+"; ";
				}
			}
			var string_final = trim(string_allData);
			var colorA = document.getElementById("homolog_a_color").value;
			var colorB = document.getElementById("homolog_b_color").value;

			var autoSubmitForm = document.createElement('form');
			autoSubmitForm.setAttribute('method','post');
			autoSubmitForm.setAttribute('action','hapmap.update_3.php');
			var input1 = document.createElement('input');
			    input1.setAttribute('type','hidden');
			    input1.setAttribute('name','hapmap_description');
			    input1.setAttribute('value',string_final);  // String containing the hapmap entries defined by the user.
			    autoSubmitForm.appendChild(input1);
			var input2 = document.createElement('input');
			    input2.setAttribute('type','hidden');
			    input2.setAttribute('name','hapmap');
			    input2.setAttribute('value','<?php echo $hapmap; ?>');
			    autoSubmitForm.appendChild(input2);
			var input3 = document.createElement('input');
			    input3.setAttribute('type','hidden');
			    input3.setAttribute('name','genome');
			    input3.setAttribute('value','<?php echo $genome; ?>');
			    autoSubmitForm.appendChild(input3);
			var input5 = document.createElement('input');
			    input5.setAttribute('type','hidden');
			    input5.setAttribute('name','project1');
			    input5.setAttribute('value','<?php echo $project1; ?>');
			    autoSubmitForm.appendChild(input5);
			var input6 = document.createElement('input');
			    input6.setAttribute('type','hidden');
			    input6.setAttribute('name','project2');
			    input6.setAttribute('value','<?php echo $project2; ?>');
			    autoSubmitForm.appendChild(input6);
			var input7 = document.createElement('input');
			    input7.setAttribute('type','hidden');
			    input7.setAttribute('name','homolog_a_color');
			    input7.setAttribute('value',colorA);
			    autoSubmitForm.appendChild(input7);
			var input8 = document.createElement('input');
			    input8.setAttribute('type','hidden');
			    input8.setAttribute('name','homolog_b_color');
			    input8.setAttribute('value',colorB);
			    autoSubmitForm.appendChild(input8);
			document.body.appendChild(autoSubmitForm);
			autoSubmitForm.submit();
		}
		function trim(str) {
			return str.replace(/^\s\s*/, '').replace(/\s\s*$/, '');
		}
	</script>
</div>
</td><td valign="top">
	<div class="tab">
		Once you've described any LOH regions for this experimental dataset, click the '<i><b>Process haplotype entry...</b></i>' button to process this entry into the haplotype map.";
	</div>
</td><tr></table>
</BODY>
</HTML>
