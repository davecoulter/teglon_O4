USE teglon;
DELIMITER $$
CREATE PROCEDURE `BackupTables`(IN map_id INT)
BEGIN
    CREATE TABLE HealpixMap_bak AS
        SELECT * FROM HealpixMap WHERE id != map_id;

    CREATE TABLE HealpixPixel_bak AS
            SELECT * FROM HealpixPixel WHERE HealpixMap_id != map_id;

    CREATE TABLE HealpixPixel_Completeness_bak AS
            SELECT * FROM HealpixPixel_Completeness WHERE HealpixMap_id != map_id;

    CREATE TABLE HealpixPixel_Galaxy_Weight_bak AS
            SELECT * FROM HealpixPixel_Galaxy_Weight WHERE HealpixMap_id != map_id;

    CREATE TABLE ObservedTile_bak AS
            SELECT * FROM ObservedTile WHERE HealpixMap_id != map_id;

    CREATE TABLE ObservedTile_HealpixPixel_bak AS
            SELECT * FROM ObservedTile_HealpixPixel WHERE HealpixMap_id != map_id;

    CREATE TABLE StaticTile_HealpixPixel_bak AS
            SELECT * FROM StaticTile_HealpixPixel WHERE HealpixMap_id != map_id;
END$$
DELIMITER ;