:: *** bin\Release フォルダのバックアップ処理を追加 ***
if exist "bin\Release" (
    echo.
    echo --- [バックアップ] bin\Release フォルダを %BACKUP_DIR% に移動します ---
    ren bin\Release bin\Release_old
    move bin\Release_old %BACKUP_DIR%
    echo バックアップが完了しました。
)