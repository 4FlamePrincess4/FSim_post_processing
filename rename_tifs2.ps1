Get-ChildItem -File | Where-Object { $_.Name -notlike "$($_.Directory.Name)_*" } | ForEach-Object {
    Rename-Item -Path $_.FullName -NewName "$($_.Directory.Name)_$($_.Name)"
}
