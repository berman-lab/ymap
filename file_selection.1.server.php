<?php
	session_start();
	$user = $_SESSION['user'];

	$key               = filter_input(INPUT_POST, "key",              FILTER_SANITIZE_STRING);
	$user              = filter_input(INPUT_POST, "user",             FILTER_SANITIZE_STRING);
	$project           = filter_input(INPUT_POST, "project",          FILTER_SANITIZE_STRING);
	$conclusion_script = filter_input(INPUT_POST, "script",           FILTER_SANITIZE_STRING);
	$file1_key         = filter_input(INPUT_POST, "File_selection_1", FILTER_SANITIZE_STRING);


$FTPdrop_dir   = "FTP_drop/".$user."/";
$FTPdrop_files = array_diff(glob($FTPdrop_dir."*"), array('..', '.'));
// Sort files by date, newest first.
array_multisort(SORT_DESC, $FTPdrop_files);
// Pull filename from list.
$file1_name    = $FTPdrop_files[$file1_key];

// Copy file from FTPdrop to: [Ymap_root]/users/[user]/projects/[project]/
$FTPdrop_file1 = str_replace($FTPdrop_dir,"",$file1_name);

copy($file1_name, "users/".$user."/projects/".$project."/".$FTPdrop_file1);



?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
	"http://www.w3.org/TR/html4/loose.dtd">
<BODY>
<script type="text/javascript">
	onload=function() {
		// Generate a form to post data to loaded page.
		var conclusion = document.createElement("form");
		conclusion.setAttribute("method","post");
		conclusion.setAttribute("action","<?php echo base64_decode($conclusion_script); ?>");
		var input1 = document.createElement("input");
		input1.setAttribute("type","hidden");
		input1.setAttribute("name","user");
		input1.setAttribute("value","<?php echo $user; ?>");
		conclusion.appendChild(input1);
		var input2 = document.createElement("input");
		input2.setAttribute("type","hidden");
		input2.setAttribute("name","project");
		input2.setAttribute("value","<?php echo $project; ?>");
		conclusion.appendChild(input2);
		var input3 = document.createElement("input");
		input3.setAttribute("type","hidden");
		input3.setAttribute("name","key");
		input3.setAttribute("value","<?php echo $key; ?>");
		conclusion.appendChild(input3);
		var input4 = document.createElement("input");
		input4.setAttribute("type","hidden");
		input4.setAttribute("name","fileName");
		input4.setAttribute("value","<?php echo $FTPdrop_file1; ?>");
		conclusion.appendChild(input4);

		// Automatically submit constructed form to post data to page.
		conclusion.submit();            // works in Safari, but not Firefox.
	}
</script>
</BODY>
</HTML>
