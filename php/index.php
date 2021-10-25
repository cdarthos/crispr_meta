<?php


if(isset($_POST['submit'])){

    $countfiles = count($_FILES['file']['name']);

    $dirtemp = $_POST['nom'];
    #echo "$dirtemp\n";
    preg_match_all('!\d+!', $dirtemp, $ok);
    #print_r($ok[0][0]);
    $dir = strval($ok[0][0]).'.'.strval($ok[0][1]).'_'.strval($ok[0][2]);

    echo $dir;
    for($i=0;$i<$countfiles;$i++){
        $filename = $_FILES['file']['name'][$i];
        if (!file_exists('upload/'.$dir)) {    mkdir('upload/'.$dir, 0777, true); }
        move_uploaded_file($_FILES['file']['tmp_name'][$i],'upload/'.$dir.'/'.$filename);
    }
}
echo "<p>Cette premiere partie permet de deposer les articles (avec les \"supplementary file\") liés au CRISPR de la Tuberculose</p>";
echo "<form method='post' action='' enctype='multipart/form-data'><p>code \"doi\" de l'article : <input type='text' name='nom'></p><input type='file' name='file[]' id='file' multiple><input type='submit' name='submit' value='Upload'></form>";

echo 'Liste des données déjà importé :<br>';
$array = array_diff(scandir('upload/'), array('..', '.'));
echo '<table>';

foreach($array as $item) {
    echo "<tr><td>$item</td>";
    $fich = array_diff(scandir('upload/'.$item), array('..', '.'));
    foreach($fich as $item1) {
        echo "<td>$item1</td>";
    }
    echo '</tr>';
}
echo '</table>';




