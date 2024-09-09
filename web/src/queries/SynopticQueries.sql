USE teglon;

# GW190425 Table 1 Column Properties
WITH UniquePixels As (
SELECT
    DISTINCT HP.id            as _2d_Pix_id,
    HPC.id           as _4d_Pix_id,
    HP.Pixel_Index,
    HP.Prob          as _2d_Prob,
    HPC.NetPixelProb as _4d_Prob,
    HP.Mean          as Mean_Dist,
    HP.Stddev        as Dist_Err,
    HPC.PixelCompleteness,
    D.Name as detect_name
FROM HealpixPixel HP
JOIN HealpixPixel_Completeness HPC on HPC.HealpixPixel_id = HP.id
JOIN ObservedTile_HealpixPixel OTHP on HP.id = OTHP.HealpixPixel_id
JOIN ObservedTile OT on OT.id = OTHP.ObservedTile_id
JOIN Detector D on OT.Detector_id = D.id
WHERE HP.HealpixMap_id = 1),
ProbSynopsis AS (
SELECT
	ut.detect_name,
	SUM(ut._2d_Prob) * 100 as sum_2d,
    SUM(ut._4d_Prob) * 100 as sum_4d
FROM 
	UniquePixels ut
GROUP BY
	ut.detect_name),
All_Tiles AS (
	SELECT
		D.`Name` as detector_name,
        D.id as detector_id,
		OT.RA,
		OT._Dec,
		OT.MJD,
		B.`Name` as filter_name,
		OT.Mag_Lim
	FROM ObservedTile OT
	JOIN Detector D on D.id = OT.Detector_id
	JOIN Band B on B.id = OT.Band_id
	WHERE OT.HealpixMap_id = 1),
GroupedTiles AS (
	SELECT 
		_at.detector_name,
		D.Area as FOV,
		COUNT(_at.MJD) AS TotalTiles
	FROM All_Tiles _at
	JOIN Detector D on D.id = _at.detector_id
	GROUP BY
		_at.detector_name, 
		D.Area
	ORDER BY
		D.Area
)
SELECT 
	CONCAT(ps.detect_name, " & ", ROUND(gt.FOV, 4), " & ", gt.TotalTiles, " & ", ROUND(ps.sum_2d, 2), " & ", ROUND(ps.sum_4d, 2), " & ", ROUND(ps.sum_4d/ps.sum_2d, 2), " \\\\") as row_tex,
    ps.detect_name,
    gt.FOV,
    gt.TotalTiles,
    ps.sum_2d, 
    ps.sum_4d, 
    ps.sum_4d/ps.sum_2d as eff_boost
FROM ProbSynopsis ps
JOIN GroupedTiles gt on ps.detect_name = gt.detector_name
ORDER BY gt.FOV;
    
    
# SQ Deg per pixel: 0.052455852825697924

WITH All_Unique AS (
SELECT
    DISTINCT HP.id            as _2d_Pix_id,
    HPC.id           as _4d_Pix_id,
    HP.Pixel_Index,
    HP.Prob          as _2d_Prob,
    HPC.NetPixelProb as _4d_Prob
FROM HealpixPixel HP
JOIN HealpixPixel_Completeness HPC on HPC.HealpixPixel_id = HP.id
JOIN ObservedTile_HealpixPixel OTHP on HP.id = OTHP.HealpixPixel_id
WHERE HP.HealpixMap_id = 1)
SELECT
	COUNT(_2d_Pix_id) * 0.052455852825697924 as TotalArea,# Area of 
	SUM(_2d_Prob)*100,
    SUM(_4d_Prob)*100
FROM All_Unique;


SELECT 
	COUNT(*)
FROM ObservedTile 
WHERE 
	HealpixMap_id = 1 AND
    RA >= 95 AND RA < 315; # 10363
    
SELECT 
	COUNT(*) 
FROM ObservedTile 
WHERE 
	HealpixMap_id = 1 AND
    RA >= 315 AND RA <= 365; # 266

SELECT 
	COUNT(*) 
FROM ObservedTile 
WHERE 
	HealpixMap_id = 1 AND
    RA > 0 AND RA <= 95; # 355
    
SELECT COUNT(*) FROM ObservedTile WHERE HealpixMap_id = 1; #10984


WITH EasternSpur AS (
SELECT 
	DISTINCT hp.id
FROM ObservedTile ot
JOIN ObservedTile_HealpixPixel ot_hp on ot_hp.ObservedTile_id = ot.id
JOIN HealpixPixel hp on hp.id = ot_hp.HealpixPixel_id
JOIN HealpixPixel_Completeness hpc on hpc.HealpixPixel_id = hp.id
WHERE 
	ot.HealpixMap_id = 1 AND
    RA >= 95 AND RA < 315)
SELECT 
	SUM(hp.Prob) as SUM2d, # 0.4516883968054497
    SUM(hpc.NetPixelProb) as SUM4d # 0.448283291570691
FROM HealpixPixel hp 
JOIN HealpixPixel_Completeness hpc on hp.id = hpc.HealpixPixel_id
JOIN EasternSpur es on es.id = hp.id;

WITH WesternSpur AS (
SELECT 
	DISTINCT hp.id
FROM ObservedTile ot
JOIN ObservedTile_HealpixPixel ot_hp on ot_hp.ObservedTile_id = ot.id
JOIN HealpixPixel hp on hp.id = ot_hp.HealpixPixel_id
JOIN HealpixPixel_Completeness hpc on hpc.HealpixPixel_id = hp.id
WHERE 
	ot.HealpixMap_id = 1 AND
    RA > 0 AND RA <= 95)
SELECT 
	SUM(hp.Prob) as SUM2d, # 0.031345386354921366
    SUM(hpc.NetPixelProb) as SUM4d # 0.03566091034261245
FROM HealpixPixel hp 
JOIN HealpixPixel_Completeness hpc on hp.id = hpc.HealpixPixel_id
JOIN WesternSpur ws on ws.id = hp.id;
    
    
WITH NotShown AS (
SELECT 
	DISTINCT hp.id
FROM ObservedTile ot
JOIN ObservedTile_HealpixPixel ot_hp on ot_hp.ObservedTile_id = ot.id
JOIN HealpixPixel hp on hp.id = ot_hp.HealpixPixel_id
JOIN HealpixPixel_Completeness hpc on hpc.HealpixPixel_id = hp.id
WHERE 
	ot.HealpixMap_id = 1 AND 
    RA >= 315 AND RA <= 365)
SELECT 
	SUM(hp.Prob) as SUM2d, # 0.0002465625670558152
    SUM(hpc.NetPixelProb) as SUM4d # 0.00019687880203980314
FROM HealpixPixel hp 
JOIN HealpixPixel_Completeness hpc on hp.id = hpc.HealpixPixel_id
JOIN NotShown ns on ns.id = hp.id;
    