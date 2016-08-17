<?php

// return the current size in GB of the user folder
function getUserUsageSize($userName)
{
	return shell_exec("find " . "users/".$userName . "/  -type f -iname 'complete.txt' | sed -e \"s/complete.txt//g\" | xargs du -scm | awk 'END{print $1}'") / (1000);
}

// return the size of the user quota in GB
function getUserQuota($userName)
{
	// load hardcoded quota from constants
	require_once 'constants.php';
	// check if user has a personal quota if so overriding quota
	if (file_exists("users/".$userName . "/quota.txt"))
		$quota = trim(file_get_contents("users/".$userName . "/quota.txt"));
	return $quota;
}

?>