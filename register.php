<?php
	session_start();
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
        "http://www.w3.org/TR/html4/loose.dtd">

<html lang="en">
	<head>
		<style type="text/css">
			body {font-family: arial;}
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>Ymap | New User Registration</title>
	</head>
	<body>
		<b>Ymap : New user registration.</b><br><br>
		The username doesn't need to correspond to your identity and passwords are encrypted, making them unavailable even to system administrators.<br>
		The true identity of the users can only be revealed if they choose to notify the administrators, which is generally not needed.<br><br>
		If a project doesn't seem to complete after a period of time (within a few hours), a comment requesting clarification can be left on the bug-tracking feature found in the "System" tab.<br>
		<br>
		<form action="php/register_server.php" method="post">
<!--		<b>Required:</b><br>  --!>
			<label for="user">Username: </label>                                        <input type="text"     id="user"                      name="user"><br>
			<label for="pw1">Password: </label>                                         <input type="password" id="pwOrig"                    name="pwOrig"><br>
			<label for="pw2">Retype Password: </label>                                  <input type="password" id="pwCopy"                    name="pwCopy"><br>
			<br>
<!--		<b>Optional:</b><br>
			<label for="primaryInvestigator_name">Primary Investigator Name: </label>   <input type="text"     id="primaryInvestigator_name"  name="primaryInvestigator_name"><br>
			<label for="primaryInvestigator_email">Primary Investigator Email: </label> <input type="text"     id="primaryInvestigator_email" name="primaryInvestigator_email"><br>
			<label for="researchInstitution">Research institution: </label>             <input type="text"     id="researchInstitution"       name="researchInstitution"><br>
			<label for="secondary_name">Secondary Name: </label>                        <input type="text"     id="secondary_name"            name="secondary_name"><br>
			<label for="secondary_email">Secondary Email: </label>                      <input type="text"     id="secondary_email"           name="secondary_email"><br>
			<br>   --!>
			<button type="submit">Register new user.</button>
		</form>
		<button type="button" onclick="window.location.href='index.php'">Back to Home</button>
	</body>
</html>
