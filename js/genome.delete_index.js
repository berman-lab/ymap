function deleteGenomeConfirmation(genome,key){
	$("#g2_"+genome+"_delete").html("<b><font color=\"red\">[Are you sure?]</font><button type='button' onclick='deleteGenome_yes(\""+genome+"\")'>Yes, delete.</button><button type='button' onclick='deleteGenome_no(\""+genome+"\",\""+key+"\")'>No, cancel</button></b>");
	document.getElementById(key).style.display = 'none';
}

// Requires javascript value user to be defined globally
function deleteGenome_yes(genome){
	$.ajax({
		url : 'php/genome.delete_server.php',
		type : 'post',
		data : {
			user: user,
			genome: genome
		},
		success : function(answer){
			if(answer == "COMPLETE"){
				window.location.href=window.location.href;
			}
		}
	});
}
function deleteGenome_no(genome,key){
	$("#g2_"+genome+"_delete").html("");
	document.getElementById(key).style.display = 'inline';
}

