USE teglon;
DELIMITER $$
CREATE PROCEDURE `DeleteMap`(IN map_id INT)
BEGIN
    SET FOREIGN_KEY_CHECKS=0;

    TRUNCATE HealpixMap;
    INSERT INTO HealpixMap SELECT * FROM HealpixMap_bak;
    DROP TABLE HealpixMap_bak;

    TRUNCATE HealpixPixel;
    INSERT INTO HealpixPixel SELECT * FROM HealpixPixel_bak;
    DROP TABLE HealpixPixel_bak;

    TRUNCATE HealpixPixel_Completeness;
    INSERT INTO HealpixPixel_Completeness SELECT * FROM HealpixPixel_Completeness_bak;
    DROP TABLE HealpixPixel_Completeness_bak;

    TRUNCATE HealpixPixel_Galaxy_Weight;
    INSERT INTO HealpixPixel_Galaxy_Weight SELECT * FROM HealpixPixel_Galaxy_Weight_bak;
    DROP TABLE HealpixPixel_Galaxy_Weight_bak;

    TRUNCATE ObservedTile;
    INSERT INTO ObservedTile SELECT * FROM ObservedTile_bak;
    DROP TABLE ObservedTile_bak;

    TRUNCATE ObservedTile_HealpixPixel;
    INSERT INTO ObservedTile_HealpixPixel SELECT * FROM ObservedTile_HealpixPixel_bak;
    DROP TABLE ObservedTile_HealpixPixel_bak;

    TRUNCATE StaticTile_HealpixPixel;
    INSERT INTO StaticTile_HealpixPixel SELECT * FROM StaticTile_HealpixPixel_bak;
    DROP TABLE StaticTile_HealpixPixel_bak;

    SET FOREIGN_KEY_CHECKS=1;
END$$
DELIMITER ;