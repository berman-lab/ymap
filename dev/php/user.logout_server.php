<?php
session_start();
require_once 'constants.php';

session_destroy();
?>
<script type="text/javascript">
    parent.update_interface_after_logout();
	parent.location.reload();
</script>
