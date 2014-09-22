<?php
	session_start();
	if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript">reload(); </script> <?php } else { $user = $_SESSION['user']; }
	require_once 'php/constants.php';
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";
?>
<HEAD>
	<style type="text/css">
		body {font-family: arial;}
	</style>
    <meta charset="UTF-8">
</HEAD>
<BODY>
<script type="text/javascript">
	onload=function() {
		// Generate a form to post data to loaded page.
		var conclusion = document.createElement("form");
		conclusion.setAttribute("method","post");
		conclusion.setAttribute("action","php/genome.working_server.php");

		var input1 = document.createElement("input");
		input1.setAttribute("type","hidden");
		input1.setAttribute("name","user");
		input1.setAttribute("value",user);
		conclusion.appendChild(input1);
		var input2 = document.createElement("input");
		input2.setAttribute("type","hidden");
		input2.setAttribute("name","genome");
		input2.setAttribute("value",genome);
		conclusion.appendChild(input2);
		var input3 = document.createElement("input");
		input3.setAttribute("type","hidden");
		input3.setAttribute("name","key");
		input3.setAttribute("value",key);
		conclusion.appendChild(input3);

		// Automatically submit constructed form to post data to page.
		conclusion.submit();		// works in Safari, but not Firefox.
	}
</script>
</BODY>
</HTML>
