/* ABIF chromatogram reader and peak-ratio heuristic; no network access required. */
(function (root, factory) {
    const api = factory();
    if (typeof module === "object" && module.exports) module.exports = api;
    else root.SangerABIF = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
    "use strict";
    const BASES = ["A", "C", "G", "T"];
    const CODES = {AG:"R", CT:"Y", CG:"S", AT:"W", GT:"K", AC:"M"};
    const MAX_BYTES = 32 * 1024 * 1024;

    function parseAbif(buffer) {
        if (!(buffer instanceof ArrayBuffer) || buffer.byteLength < 34 ||
                buffer.byteLength > MAX_BYTES) {
            throw new Error("Expected an AB1 file between 34 bytes and 32 MiB");
        }
        const view = new DataView(buffer);
        function check(offset, count) {
            if (!Number.isSafeInteger(offset) || offset < 0 ||
                    !Number.isSafeInteger(count) || count < 0 ||
                    offset > buffer.byteLength - count) {
                throw new Error("ABIF directory or data extends beyond the file");
            }
        }
        function str(offset, count) {
            check(offset, count);
            let result = "";
            for (let i = 0; i < count; i++) result += String.fromCharCode(view.getUint8(offset + i));
            return result;
        }
        if (str(0, 4) !== "ABIF") throw new Error("Not an ABIF chromatogram");
        function entry(offset) {
            check(offset, 28);
            return {
                name: str(offset, 4),
                num: view.getUint32(offset + 4, false),
                type: view.getUint16(offset + 8, false),
                size: view.getUint16(offset + 10, false),
                count: view.getUint32(offset + 12, false),
                bytes: view.getUint32(offset + 16, false),
                location: view.getUint32(offset + 20, false),
                offset: offset
            };
        }
        const directory = entry(6);
        if (directory.count === 0 || directory.count > 10000) throw new Error("Invalid ABIF tag count");
        check(directory.location, directory.count * 28);
        const tags = new Map();
        for (let i = 0; i < directory.count; i++) {
            const item = entry(directory.location + i * 28);
            tags.set(item.name + item.num, item);
        }
        function data(tag, fallback) {
            const item = tags.get(tag) || (fallback && tags.get(fallback));
            if (!item) return null;
            if (!item.size || item.count > MAX_BYTES ||
                    item.bytes < item.count * item.size) {
                throw new Error("Invalid ABIF tag size: " + tag);
            }
            const start = item.bytes <= 4 ? item.offset + 20 : item.location;
            check(start, item.bytes);
            return {item, start};
        }
        function chars(tag, fallback) {
            const d = data(tag, fallback);
            return d ? str(d.start, d.item.count).replace(/\0+$/, "") : null;
        }
        function numbers(tag, fallback) {
            const d = data(tag, fallback);
            if (!d) return null;
            const {item, start} = d;
            if (![1, 2, 4].includes(item.size)) throw new Error("Unexpected sample width for " + tag);
            const values = new Array(item.count);
            for (let i = 0; i < item.count; i++) {
                const offset = start + i * item.size;
                if (item.size === 1) values[i] = view.getUint8(offset);
                else if (item.size === 2) values[i] = item.type === 4 ?
                    view.getInt16(offset, false) : view.getUint16(offset, false);
                else values[i] = item.type === 5 ?
                    view.getInt32(offset, false) : view.getUint32(offset, false);
            }
            return values;
        }
        const order = chars("FWO_1");
        if (!order || order.length !== 4 || new Set(order.split("")).size !== 4 ||
                !BASES.every(b => order.includes(b))) {
            throw new Error("Missing or invalid four-dye order (FWO_1)");
        }
        const channels = {};
        for (let i = 0; i < 4; i++) {
            const values = numbers("DATA" + (i + 9));
            if (!values || values.length === 0) throw new Error("Missing analyzed trace DATA" + (i + 9));
            channels[order[i]] = values;
        }
        const length = channels.A.length;
        if (!BASES.every(b => channels[b].length === length)) {
            throw new Error("Analyzed trace channels have inconsistent lengths");
        }
        const positions = numbers("PLOC2", "PLOC1");
        const sequence = chars("PBAS2", "PBAS1");
        const qualities = numbers("PCON2", "PCON1");
        if (!positions || !sequence || positions.length !== sequence.length ||
                (qualities && qualities.length !== positions.length)) {
            throw new Error("Called bases, quality scores and peak locations do not agree");
        }
        if (!positions.every(p => p >= 0 && p < length)) {
            throw new Error("Peak position outside trace bounds");
        }
        if (!/^[ACGTN]+$/i.test(sequence)) {
            throw new Error("Unsupported base calls in PBAS");
        }
        return {sequence: sequence.toUpperCase(), positions, qualities: qualities || positions.map(() => 0), channels};
    }

    function callPeaks(trace, threshold, window, minQuality) {
        if (!Number.isFinite(threshold) || threshold <= 0 || threshold > 1) {
            throw new Error("Peak ratio threshold must be between 0 and 1");
        }
        if (!Number.isInteger(window) || window < 0 || window > 50) {
            throw new Error("Peak window must be an integer between 0 and 50");
        }
        if (!Number.isInteger(minQuality) || minQuality < 0 || minQuality > 100) {
            throw new Error("Minimum quality must be between 0 and 100");
        }
        const rows = [];
        for (let i = 0; i < trace.positions.length; i++) {
            if (trace.qualities[i] < minQuality) continue;
            const pos = trace.positions[i];
            const lo = Math.max(0, pos - window);
            const hi = Math.min(trace.channels.A.length, pos + window + 1);
            const heights = BASES.map(base => ({
                base, height: Math.max(0, ...trace.channels[base].slice(lo, hi))
            })).sort((a, b) => b.height - a.height);
            const primary = heights[0];
            const secondary = heights[1];
            const ratio = primary.height > 0 ? secondary.height / primary.height : 0;
            const flagged = secondary.height > 0 && ratio >= threshold;
            rows.push({
                position: i + 1, traceIndex: pos, calledBase: trace.sequence[i],
                quality: trace.qualities[i], primaryBase: primary.base, primaryHeight: primary.height,
                secondaryBase: secondary.base, secondaryHeight: secondary.height, ratio,
                iupac: flagged ? CODES[[primary.base, secondary.base].sort().join("")] : primary.base,
                flagged
            });
        }
        return rows;
    }

    function demoTrace() {
        const sequence = "ACGTACGTACGT";
        const positions = [...sequence].map((_, i) => 12 + 12 * i);
        const channels = Object.fromEntries(BASES.map(b => [b, Array(170).fill(15)]));
        positions.forEach((p, i) => { channels[sequence[i]][p] = 500; });
        channels.G[positions[4]] = 340;
        channels.T[positions[8]] = 180;
        return {sequence, positions, qualities: positions.map(() => 40), channels};
    }
    return {parseAbif, callPeaks, demoTrace};
});
